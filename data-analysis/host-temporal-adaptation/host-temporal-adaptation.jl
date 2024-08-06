#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
# using BenchmarkTools


run_id = ARGS[1]
delayInc = parse(Int64, ARGS[2])
sampleInc = parse(Int64, ARGS[3])
if delayInc < sampleInc
    delayInc = sampleInc
end

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("host-temporal-adaptation_output.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/26_MOI3", "sweep_db_gathered.sqlite")
# dbOutputPath = joinpath("/Volumes/Yadgah", "host-temporal-adaptation.sqlite")
rm(dbOutputPath, force=true)

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
maxT = maximum([t for (t,) in execute(dbSim, "SELECT t_final FROM param_combos")])

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
tmpPath = joinpath(ENV["SLURM_TMPDIR"], "host-adapt-$(run_id).sqlite")
# tmpPath = joinpath("/Volumes/Yadgah/test.sqlite") # local
println("this is the temp path: $(tmpPath)")
rm(tmpPath, force=true)
dbTempSim = SQLite.DB(tmpPath)
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")
execute(dbTempSim, "CREATE TABLE vpspacers (vstrain_id INTEGER, spacer_id INTEGER)")
execute(dbTempSim, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim, "INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTempSim, "INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTempSim, "INSERT INTO bspacers (bstrain_id, spacer_id) SELECT bstrain_id, spacer_id FROM dbSim.bspacers WHERE run_id = $(run_id);")
execute(dbTempSim, "INSERT INTO vpspacers (vstrain_id, spacer_id) SELECT vstrain_id, spacer_id FROM dbSim.vpspacers WHERE run_id = $(run_id);")
execute(dbTempSim, "INSERT INTO summary (t, microbial_abundance, viral_abundance) SELECT t, microbial_abundance, viral_abundance FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTempSim, "COMMIT")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id)")
execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTempSim, "CREATE INDEX bspacers_index ON bspacers (bstrain_id)")
execute(dbTempSim, "CREATE INDEX vpspacers_index ON vpspacers (vstrain_id)")
execute(dbTempSim, "CREATE INDEX summary_index ON summary (t)")
execute(dbTempSim, "COMMIT")


function virusesWithMissingTimes(tDelays, sampleInc)
    vtimes = [t for (t,) in execute(dbTempSim, "SELECT t FROM summary")]
    tMextinct = [t for (t,) in execute(dbTempSim, "SELECT t
                FROM summary 
                WHERE viral_abundance == 0")]
    if length(tMextinct) > 0 
        if minimum(tMextinct) < maximum(tDelays)
            tMax = maximum(vtimes[vtimes.<minimum(tMextinct)])
            times = collect(0:sampleInc:tMax)
        else
            times = tDelays
        end
    else
        times = tDelays
    end
    vstrains = DataFrame(execute(
            dbTempSim,
            "SELECT DISTINCT t, viral_abundance FROM summary 
                WHERE t in ($(join(times,", ")))"))
    vstrains = innerjoin(vstrains, DataFrame(execute(
                dbTempSim, "SELECT DISTINCT t, vstrain_id, abundance FROM vabundance 
                WHERE t in ($(join(times,", "))) ORDER BY vstrain_id")), on=:t)                
    timesMissing = sort([setdiff(Set(times), Set(vstrains[:,:t]))...])
    # if length(tMextinct) > 0 
    #     if tMextinct < maximum(tDelays)
    #         timesMissing = timesMissing[timesMissing .< tMextinct]
    #     end
    # end
    if length(timesMissing) > 0
        for tMiss in timesMissing
            lastTime = maximum(vtimes[vtimes.<tMiss])
            # println("last time: $(lastTime)")
            new = DataFrame(execute(
                dbTempSim,
                "SELECT DISTINCT t, viral_abundance FROM summary 
                    WHERE t = $(lastTime)"))
            new = innerjoin(new, DataFrame(execute(
                dbTempSim, "SELECT DISTINCT t, vstrain_id, abundance FROM vabundance 
                WHERE t = $(lastTime) ORDER BY vstrain_id")),on=:t)
            replace!(new.t, lastTime => tMiss)
            append!(vstrains, new)
        end
    end
    rename!(vstrains, :viral_abundance => :vfrequency)
    rename!(vstrains, :t => :t_delay)
    vstrains[!, :vfrequency] = vstrains[!, :abundance] ./ vstrains[!, :vfrequency]
    select!(vstrains, Not([:abundance]))
    return vstrains, tMextinct
end
function identifyMatches!(bstrainsDF, sampleInc, dbOutput)
    tDelays = unique(bstrainsDF[:, :t_delay])
    vstrains, tMextinct = virusesWithMissingTimes(tDelays, sampleInc)
    # if length(tMextinct) > 0
    #     if maximum(tDelays) > minimum(tMextinct)
    #         vdelays = Vector{Float64}(bstrainsDF[bstrainsDF.t_delay .>=minimum(tMextinct), :t_delay])
    #         numTimes = length(vdelays)
    #         svstrains = DataFrame(t_delay=vdelays, vfrequency=Vector{Float64}(repeat([0], numTimes)))
    #         innerjoin(bstrainsDF,svstrains,on=:t_delay) |> 
    #             SQLite.load!(dbOutput, "btemporal_adaptation", ifnotexists=true)
    #     end
    # end
    bstrainID = bstrainsDF[:,:bstrain_id][1]
    svstrains = [spacerID for (spacerID,) in execute(
        dbTempSim,
        "SELECT spacer_id FROM bspacers
        WHERE bstrain_id = $(bstrainID)")]
    strains = vstrains[:, :vstrain_id]
    svstrains = DataFrame(execute(dbTempSim, "SELECT vstrain_id, spacer_id FROM vpspacers 
        WHERE vstrain_id in ($(join(strains,", "))) AND spacer_id in ($(join(svstrains,", ")))"
    ))
    svstrains = svstrains[:, :vstrain_id]

    DataFrame(bstrain_id=repeat([bstrainID],length(svstrains)),vstrain_id=svstrains) |> 
            SQLite.load!(dbOutput, "bstrain_to_vstrain_matches", ifnotexists=true)
    svstrains = vstrains[[in(x, svstrains) for x in vstrains.vstrain_id], :]
    svstrains = groupby(svstrains, :t_delay)
    svstrains = combine(svstrains, [:vfrequency] .=> sum; renamecols=false)
    innerjoin(bstrainsDF, svstrains, on=:t_delay) |> SQLite.load!(dbOutput, "btemporal_adaptation", ifnotexists=true)
end

function temporalAdaptation(maxT, sampleInc, delayInc, dbTempSim, dbOutput)
    delays = unique([reverse(collect(0:-delayInc:-maxT))..., collect(0:delayInc:maxT)...])
    times = collect(0:sampleInc:maxT)
    map(t->repeat([t],length(delays)),times)
    timeDelayDF = DataFrame(t = reduce(vcat, map(t -> repeat([t], length(delays)), times)), delay = reduce(vcat, repeat(delays, length(times))))
    bstrains = [bstrainID for (bstrainID,) in execute(
        dbTempSim, "SELECT DISTINCT bstrain_id FROM babundance 
         WHERE t in ($(join(times,", "))) ORDER BY bstrain_id")]
    bTotalDF = DataFrame(execute(
        dbTempSim, "SELECT DISTINCT t, microbial_abundance FROM summary 
         WHERE t in ($(join(times,", ")))"))
    rename!(bTotalDF, :microbial_abundance => :bfrequency)
    for bstrainID in bstrains
        println("bstrainID: $(bstrainID)")
        bstrainsDF = innerjoin(DataFrame(execute(dbTempSim, "SELECT DISTINCT t, bstrain_id, abundance FROM babundance 
                        WHERE t in ($(join(times,", "))) AND bstrain_id = $(bstrainID)")),timeDelayDF,on=:t)
        bstrainsDF = innerjoin(bstrainsDF,bTotalDF,on=:t)
        bstrainsDF[!, :bfrequency] = bstrainsDF[!, :abundance] ./ bstrainsDF[!, :bfrequency]
        select!(bstrainsDF,Not([:abundance]))
        bstrainsDF[!, :t_delay] = bstrainsDF[!, :t] .+ bstrainsDF[!, :delay]
        bstrainsDF = bstrainsDF[(bstrainsDF.t_delay.<=maxT).&(bstrainsDF.t_delay.>=0), :]
        identifyMatches!(bstrainsDF, sampleInc, dbOutput)
    end
end

temporalAdaptation(maxT, sampleInc, delayInc, dbTempSim, dbOutput)
rm(tmpPath, force=true)
println("Complete!")
