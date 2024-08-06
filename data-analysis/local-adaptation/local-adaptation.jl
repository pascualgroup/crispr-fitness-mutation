#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = parse(Int64, ARGS[1])
# run_id2 = ARGS[2]
timeInc = parse(Int64, ARGS[2])

if length(ARGS) > 2
    threshold = parse(Float64, ARGS[3])/100
else
    threshold = 0
end
## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("local-adaptation_output.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/21_MOI3", "sweep_db_gathered.sqlite")
# dbOutputPath = joinpath("/Volumes/Yadgah/local-adaptation.sqlite") # local
rm(dbOutputPath, force=true)


# dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
# execute(dbOutput, "CREATE TABLE vlocal_adaptation (vrun_id INTEGER, brun_id INTEGER, t REAL,  vstrain_id INTEGER, vfrequency REAL, bfrequency REAL)")

tmpPath = joinpath(ENV["SLURM_TMPDIR"], "local-adapt-$(run_id).sqlite")
# tmpPath = dbSimPath
# tmpPath = joinpath("/Volumes/Yadgah/test.sqlite")
println("this is the temp path: $(tmpPath)")
dbTempSim = SQLite.DB(tmpPath)




execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER, run_id INTEGER)")
execute(dbTempSim, "CREATE TABLE vpspacers (vstrain_id INTEGER, spacer_id INTEGER, run_id INTEGER)")
execute(dbTempSim, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER, run_id INTEGER)")
execute(dbTempSim, "BEGIN TRANSACTION")
execute(dbTempSim, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTempSim, "INSERT INTO vabundance (t, vstrain_id, abundance, run_id) SELECT t, vstrain_id, abundance, run_id FROM dbSim.vabundance WHERE run_id in ($(run_id));")
execute(dbTempSim, "INSERT INTO vpspacers (vstrain_id, spacer_id, run_id) SELECT vstrain_id, spacer_id, run_id FROM dbSim.vpspacers WHERE run_id in ($(run_id));")
execute(dbTempSim, "INSERT INTO summary (t, microbial_abundance, viral_abundance, run_id) SELECT t, microbial_abundance, viral_abundance, run_id FROM dbSim.summary WHERE run_id in ($(run_id));")
execute(dbTempSim, "COMMIT")
execute(dbTempSim, "DETACH DATABASE dbSim")

execute(dbTempSim, "BEGIN TRANSACTION")

execute(dbTempSim, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTempSim, "CREATE INDEX vpspacers_index ON vpspacers (vstrain_id)")
execute(dbTempSim, "CREATE INDEX summary_index ON summary (t)")
execute(dbTempSim, "COMMIT")

# execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER, run_id INTEGER)")
# execute(dbTempSim, "INSERT INTO babundance (t, bstrain_id, abundance, run_id) SELECT t, bstrain_id, abundance, run_id FROM dbSim.babundance WHERE run_id in ($(run_id));")
# execute(dbTempSim, "CREATE INDEX babundance_index ON babundance (t,bstrain_id)")
# execute(dbTempSim, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER, run_id INTEGER)")
# execute(dbTempSim, "INSERT INTO bspacers (bstrain_id, spacer_id, run_id) SELECT bstrain_id, spacer_id, run_id FROM dbSim.bspacers WHERE run_id in ($(run_id));")
# execute(dbTempSim, "CREATE INDEX bspacers_index ON bspacers (bstrain_id)")

dbSim = SQLite.DB(dbSimPath)
(cID,) = execute(dbSim, "SELECT combo_id FROM runs WHERE run_id = $(run_id)")
cID = cID.combo_id
runIDs2 = [runID for (runID,) in execute(dbSim, "SELECT run_id FROM runs WHERE combo_id = $(cID)")]


function microbesWithMissingTimes(dbSim, run_id2, vtimes, timeInc)
    btimes = [t for (t,) in 
            execute(dbSim, "SELECT t FROM summary 
            WHERE run_id = $(run_id2)")]
    tMextinct = [t for (t,) in execute(dbSim, "SELECT t
                FROM summary 
                WHERE run_id = $(run_id2) AND microbial_abundance == 0")]
    if length(tMextinct) > 0 
        if tMextinct < maximum(vtimes)
            tMax = maximum(btimes[btimes.<tMextinct])
            times = collect(0:timeInc:tMax)
        else
            times = copy(vtimes[:, :t])
        end
    else
        times = copy(vtimes[:, :t])
    end  
    bfreq = DataFrame(execute(dbSim, "SELECT t, microbial_abundance 
                FROM summary 
                WHERE run_id = $(run_id2)
                AND t in ($(join(times,", ")))"))
    bfreq = innerjoin(bfreq, DataFrame(execute(dbSim, "SELECT t, bstrain_id, abundance 
                FROM babundance 
                WHERE run_id = $(run_id2)
                AND t in ($(join(times,", ")))")),on=:t)
    timesMissing = sort([setdiff(Set(times),Set(bfreq[!,:t]))...])
    # if length(tMextinct) > 0
    #     if tMextinct < maximum(tDelays)
    #         timesMissing = timesMissing[timesMissing.<tMextinct]
    #     end
    # end
    if length(timesMissing) > 0
        for tMiss in timesMissing
            lastTime = maximum(btimes[btimes .< tMiss])
            new = DataFrame(execute(dbSim, "SELECT t, microbial_abundance 
                    FROM summary 
                    WHERE run_id = $(run_id2)
                    AND t = $(lastTime)"))
            new = innerjoin(new, DataFrame(execute(dbSim, "SELECT t, bstrain_id, abundance 
                    FROM babundance 
                    WHERE run_id = $(run_id2)
                    AND t = $(lastTime)")),on=:t)
            replace!(new.t,lastTime => tMiss)
            append!(bfreq,new)
        end
    end
    rename!(bfreq, :abundance => :bfrequency)
    bfreq[!, :bfrequency] = bfreq[!, :bfrequency] ./ bfreq[!, :microbial_abundance]
    select!(bfreq, Not([:microbial_abundance]))
    return bfreq, tMextinct
end


function identifyAll0Matches!(run_id, runIDs2, vstrainID, vfreq, timeInc)
    vstrains = DataFrame(execute(dbTempSim, "SELECT vstrain_id, spacer_id 
            FROM vpspacers WHERE run_id = $(run_id) 
            AND vstrain_id = $(vstrainID)"))
    for run_id2 in runIDs2
        sbstrains = copy(vstrains)
        # println("runID: $(run_id2)")
        bfreq, tMextinct = microbesWithMissingTimes(dbSim, run_id2, vfreq, timeInc)
        if length(tMextinct) > 0
            # this is if the microbe population goes extinctin. 
            # tMissing are the times the microbe is extinct but not the virus
            if maximum(vfreq[:, :t]) > tMextinct[1]
                vtimes = Vector{Float64}(vfreq[vfreq.t.>=tMextinct[1], :t])
                numTimes = length(vtimes)
                DataFrame(vrun_id=Vector{Int64}(repeat([run_id], numTimes)),
                    brun_id=Vector{Int64}(repeat([run_id2], numTimes)), t=vtimes,
                    vstrain_id=Vector{Int64}(repeat([vstrainID], numTimes)),
                    vfrequency=Vector{Float64}(repeat([vfreq], numTimes)),
                    bfrequency=Vector{Float64}(repeat([0], numTimes))) |>
                SQLite.load!(dbOutput, "vlocal_adaptation", ifnotexists=true)
            end
        end
        strains = unique(bfreq[:,:bstrain_id])
        sbstrains = DataFrame(execute(dbSim, "SELECT bstrain_id, spacer_id 
                    FROM bspacers WHERE run_id = $(run_id2) 
                    AND bstrain_id in ($(join(strains,", ")))
                    AND spacer_id in ($(join(sbstrains[:,:spacer_id],", ")))"))
        sbstrains = setdiff(strains, sbstrains[:, :bstrain_id])
        # if length(sbstrains) == 0
        #     println(vstrainID)
        #     # DataFrame(vrun_id=Vector{Int64}(repeat([run_id], numTimes)),
        #     #     brun_id=Vector{Int64}(repeat([run_id2], numTimes)), t=vtimes,
        #     #     vstrain_id=Vector{Int64}(repeat([vstrainID], numTimes)),
        #     #     vfrequency=Vector{Float64}(vfreq[:, :vfrequency]),
        #     #     bfrequency=Vector{Float64}(repeat([0], numTimes))) |>
        #     # SQLite.load!(dbOutput, "vlocal_adaptation", ifnotexists=true)
        #     return
        # end
        bfreq = bfreq[[in(x, sbstrains) for x in bfreq.bstrain_id], :]
        bfreq = groupby(bfreq, :t)
        bfreq = combine(bfreq, [:bfrequency] .=> sum; renamecols=false)
        vtimes = Vector{Float64}(vfreq[:, :t])
        numTimes = length(vtimes)
        if length(sbstrains) > 0
            innerjoin(DataFrame(vrun_id=Vector{Int64}(repeat([run_id], numTimes)),
                    brun_id=Vector{Int64}(repeat([run_id2], numTimes)), t=vtimes,
                    vstrain_id=Vector{Int64}(repeat([vstrainID], numTimes)),
                    vfrequency=Vector{Float64}(vfreq[:, :vfrequency])), 
                    bfreq, on =:t) |>
                        SQLite.load!(dbOutput, "vlocal_adaptation", ifnotexists=true)
        end
        # this is if the microbe population goes extinctin. 
        # tMissing are the times the microbe is extinct but not the virus
        # tMissing = Vector{Float64}(setdiff(vtimes, bfreq[:, :t]))
        # numTimes = length(tMissing)
        # if numTimes > 0
        #     println("Immune types exist for strain: $(vstrainID)")
        #     DataFrame(vrun_id=Vector{Int64}(repeat([run_id], numTimes)),
        #             brun_id=Vector{Int64}(repeat([run_id2], numTimes)), t=tMissing,
        #             vstrain_id=Vector{Int64}(repeat([vstrainID], numTimes)),
        #             vfrequency=Vector{Float64}(vfreq[[in(x, tMissing) for x in vfreq.t],:][:, :vfrequency]),
        #             bfrequency=Vector{Float64}(repeat([0], numTimes))) |>
        #     SQLite.load!(dbOutput, "vlocal_adaptation", ifnotexists=true)
        # end
    end
    return
end


function localAdaptation(run_id, runIDs2, timeInc, threshold)
    maxT = maximum([t for (t,) in 
            execute(dbTempSim, "SELECT t FROM summary 
            WHERE run_id = $(run_id) 
            AND viral_abundance != 0")])
    times = collect(0:timeInc:maxT)
    vfreq = DataFrame(execute(dbTempSim, "SELECT t, viral_abundance 
            FROM summary WHERE run_id = $(run_id)
            AND t in ($(join(times,", ")))"
    ))
    rename!(vfreq, :viral_abundance => :vtotal)
    vfreq = innerjoin(vfreq,DataFrame(execute(dbTempSim, "SELECT t, vstrain_id, abundance 
            FROM vabundance WHERE run_id = $(run_id)
            AND t in ($(join(times,", ")))"
        )), on=:t)
    rename!(vfreq, :abundance => :vfrequency)
    vfreq[!,:vfrequency] =  vfreq[!,:vfrequency]./vfreq[!,:vtotal]
    select!(vfreq, Not([:vtotal]))
    filteredStrains = unique(vfreq[vfreq.vfrequency .>=threshold, :][:,:vstrain_id])
    # println("checking $(run_id2)")
    #select all times here and the filter out    
  
    for vstrainID in filteredStrains
        println("vstrainID: $(vstrainID)")
        identifyAll0Matches!(run_id, runIDs2, vstrainID, vfreq[vfreq.vstrain_id.==vstrainID, :], timeInc)
    end 
end


localAdaptation(run_id,runIDs2,timeInc, threshold)
println("Complete!")