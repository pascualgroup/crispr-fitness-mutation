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
dbOutputPath = joinpath("temporal-adaptation_output.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/21_MOI3", "sweep_db_gathered.sqlite")
# dbOutputPath = joinpath("/Volumes/Yadgah", "temporal-adaptation2.sqlite")
rm(dbOutputPath, force=true)

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
maxT = maximum([t for (t,) in execute(dbSim, "SELECT t_final FROM param_combos")])

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
tmpPath = joinpath(ENV["SLURM_TMPDIR"], "temp-adapt2-$(run_id).sqlite")
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


function microbesWithMissingTimes(tDelays, sampleInc)
    btimes = [t for (t,) in execute(dbTempSim, "SELECT t FROM summary")]
    tMextinct = [t for (t,) in execute(dbTempSim, "SELECT t
                FROM summary 
                WHERE microbial_abundance == 0")]
    if length(tMextinct) > 0 
        if minimum(tMextinct) < maximum(tDelays)
            tMax = maximum(btimes[btimes.<minimum(tMextinct)])
            times = collect(0:sampleInc:tMax)
        else
            times = tDelays
        end
    else
        times = tDelays
    end
    bstrains = DataFrame(execute(
            dbTempSim,
            "SELECT DISTINCT t, microbial_abundance FROM summary 
                WHERE t in ($(join(times,", ")))"))
    bstrains = innerjoin(bstrains, DataFrame(execute(
                dbTempSim, "SELECT DISTINCT t, bstrain_id, abundance FROM babundance 
                WHERE t in ($(join(times,", "))) ORDER BY bstrain_id")), on=:t)                
    timesMissing = sort([setdiff(Set(times), Set(bstrains[:,:t]))...])
    # if length(tMextinct) > 0 
    #     if tMextinct < maximum(tDelays)
    #         timesMissing = timesMissing[timesMissing .< tMextinct]
    #     end
    # end
    if length(timesMissing) > 0
        for tMiss in timesMissing
            lastTime = maximum(btimes[btimes.<tMiss])
            # println("last time: $(lastTime)")
            new = DataFrame(execute(
                dbTempSim,
                "SELECT DISTINCT t, microbial_abundance FROM summary 
                    WHERE t = $(lastTime)"))
            new = innerjoin(new, DataFrame(execute(
                dbTempSim, "SELECT DISTINCT t, bstrain_id, abundance FROM babundance 
                WHERE t = $(lastTime) ORDER BY bstrain_id")),on=:t)
            replace!(new.t, lastTime => tMiss)
            append!(bstrains, new)
        end
    end
    rename!(bstrains, :microbial_abundance => :bfrequency)
    rename!(bstrains, :t => :t_delay)
    bstrains[!, :bfrequency] = bstrains[!, :abundance] ./ bstrains[!, :bfrequency]
    select!(bstrains, Not([:abundance]))
    return bstrains, tMextinct
end
function identify0Matches!(vstrainsDF, sampleInc, dbOutput)
    tDelays = unique(vstrainsDF[:, :t_delay])
    bstrains, tMextinct = microbesWithMissingTimes(tDelays, sampleInc)
    if length(tMextinct) > 0
        if maximum(tDelays) > minimum(tMextinct)
            vdelays = Vector{Float64}(vstrainsDF[vstrainsDF.t_delay .>=minimum(tMextinct), :t_delay])
            numTimes = length(vdelays)
            sbstrains = DataFrame(t_delay=vdelays, bfrequency=Vector{Float64}(repeat([0], numTimes)))
            innerjoin(vstrainsDF,sbstrains,on=:t_delay) |> 
                SQLite.load!(dbOutput, "vtemporal_adaptation", ifnotexists=true)
        end
    end
    vstrainID = vstrainsDF[:,:vstrain_id][1]
    sbstrains = [spacerID for (spacerID,) in execute(
        dbTempSim,
        "SELECT spacer_id FROM vpspacers
        WHERE vstrain_id = $(vstrainID)")]
    strains = bstrains[:,:bstrain_id]
    sbstrains = DataFrame(execute(dbTempSim, "SELECT bstrain_id, spacer_id FROM bspacers 
        WHERE bstrain_id in ($(join(strains,", "))) AND spacer_id in ($(join(sbstrains,", ")))"
    ))
    sbstrains = setdiff(strains, sbstrains[:, :bstrain_id])
    # if length(sbstrains) > 0
    DataFrame(vstrain_id=repeat([vstrainID],length(sbstrains)),bstrain_id=sbstrains) |> 
            SQLite.load!(dbOutput, "bstrain_to_vstrain_0matches", ifnotexists=true)
    sbstrains = bstrains[[in(x, sbstrains) for x in bstrains.bstrain_id], :]
    sbstrains = groupby(sbstrains, :t_delay)
    sbstrains = combine(sbstrains, [:bfrequency] .=> sum; renamecols=false)
    innerjoin(vstrainsDF, sbstrains, on=:t_delay) |> SQLite.load!(dbOutput, "vtemporal_adaptation", ifnotexists=true)
    # end
    # tMissing = Vector{Float64}(setdiff(tDelays, sbstrains[:, :t_delay]))
    # numTimes = length(tMissing)
    # if numTimes > 0
    #     # println("Immune types exist for strain: $(vstrainID)")
    #     innerjoin(vstrainsDF[[in(x, tMissing) for x in vstrainsDF.t_delay],:], 
    #     DataFrame(t_delay=tMissing, bfrequency=Vector{Float64}(repeat([0],numTimes))), 
    #         on=:t_delay) |> SQLite.load!(dbOutput, "vtemporal_adaptation", ifnotexists=true)
    # end
end

function temporalAdaptation(maxT, sampleInc, delayInc, dbTempSim, dbOutput)
    delays = unique([reverse(collect(0:-delayInc:-maxT))..., collect(0:delayInc:maxT)...])
    times = collect(0:sampleInc:maxT)
    map(t->repeat([t],length(delays)),times)
    timeDelayDF = DataFrame(t = reduce(vcat, map(t -> repeat([t], length(delays)), times)), delay = reduce(vcat, repeat(delays, length(times))))
    vstrains = [vstrainID for (vstrainID,) in execute(
        dbTempSim, "SELECT DISTINCT vstrain_id FROM vabundance 
         WHERE t in ($(join(times,", "))) ORDER BY vstrain_id")]
    vTotalDF = DataFrame(execute(
        dbTempSim, "SELECT DISTINCT t, viral_abundance FROM summary 
         WHERE t in ($(join(times,", ")))"))
    rename!(vTotalDF, :viral_abundance => :vfrequency)
    for vstrainID in vstrains
        println("vstrainID: $(vstrainID)")
        vstrainsDF = innerjoin(DataFrame(execute(dbTempSim, "SELECT DISTINCT t, vstrain_id, abundance FROM vabundance 
                        WHERE t in ($(join(times,", "))) AND vstrain_id = $(vstrainID)")),timeDelayDF,on=:t)
        vstrainsDF = innerjoin(vstrainsDF,vTotalDF,on=:t)
        vstrainsDF[!, :vfrequency] = vstrainsDF[!, :abundance] ./ vstrainsDF[!, :vfrequency]
        select!(vstrainsDF,Not([:abundance]))
        vstrainsDF[!, :t_delay] = vstrainsDF[!, :t] .+ vstrainsDF[!, :delay]
        vstrainsDF = vstrainsDF[(vstrainsDF.t_delay.<=maxT).&(vstrainsDF.t_delay.>=0), :]
        identify0Matches!(vstrainsDF, sampleInc, dbOutput)
    end
end

temporalAdaptation(maxT, sampleInc, delayInc, dbTempSim, dbOutput)
rm(tmpPath, force=true)
println("Complete!")
