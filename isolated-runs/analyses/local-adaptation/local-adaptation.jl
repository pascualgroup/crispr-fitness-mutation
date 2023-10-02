#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
# using DataFramesMeta
using SQLite.DBInterface: execute
# using BenchmarkTools


run_id = ARGS[1]
timeInc = parse(Int64, ARGS[2])

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("local-adaptation_output.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/1_MOI3-vary-fgm-scale/simulation", "sweep_db_gathered.sqlite")
# dbTriPath = joinpath("/Volumes/Yadgah/comboID66/tripartite-networksC66.sqlite") # local
# dbShanPath = joinpath("/Volumes/Yadgah/comboID66/shannonC66.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah/comboID66/matchesC66.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/comboID66/temporal-adaptation.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/comboID1/local-adaptation.sqlite") # local

dbTempSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
cID = [id for (id,) in execute(dbTempSim, "SELECT combo_id FROM 
        runs WHERE run_id = $(run_id)")]
run_ids2 = [run_id for (run_id,) in 
            execute(dbTempSim, "SELECT run_id FROM 
            runs WHERE combo_id = $(cID[1])")]

execute(dbOutput, "CREATE TABLE local_adaptation (brun_id INTEGER, t REAL,  vstrain_id INTEGER, vfitness REAL)")

function localAdaptation(run_id,run_ids2,timeInc)
    vstrains = [ID for (ID,) in execute(dbTempSim, 
                "SELECT DISTINCT vstrain_id 
                FROM vabundance WHERE run_id = $(run_id)")]
    vfreq = DataFrame(execute(dbTempSim, "SELECT t, viral_abundance 
            FROM summary WHERE run_id = $(run_id)
            AND viral_abundance != 0"))
    vfreq = vfreq[vfreq.t .% timeInc .== 0 ,:]
    rename!(vfreq, :viral_abundance => :vtotal)
    vfreq = innerjoin(vfreq,DataFrame(execute(dbTempSim, "SELECT t, vstrain_id, abundance 
            FROM vabundance WHERE run_id = $(run_id)
            AND t in ($(join(unique(vfreq[!,:t]),", ")))")),on=:t)
    rename!(vfreq, :abundance => :vfreq)
    vfreq[!,:vfreq] =  vfreq[!,:vfreq]./vfreq[!,:vtotal]
    select!(vfreq, Not([:vtotal]))
    vspacers = innerjoin(vfreq[!,[:vstrain_id]],DataFrame(execute(dbTempSim, "SELECT vstrain_id, spacer_id 
            FROM vpspacers WHERE run_id = $(run_id)")),on=:vstrain_id)
    for (run_id2,) in run_ids2
        println("checking $(run_id2)")
        (maxT,) = execute(dbTempSim, "SELECT DISTINCT t_final FROM param_combos") 
        maxT = maxT.t_final
        maxT = maxT - (maxT % timeInc)
        times = collect(0:timeInc:maxT)
        btotal = DataFrame(execute(dbTempSim, "SELECT t, microbial_abundance 
                FROM summary 
                WHERE run_id = $(run_id2) AND microbial_abundance != 0
                AND t in ($(join(times,", ")))"))
        rename!(btotal, :microbial_abundance => :btotal)
        timesMissing = sort([setdiff(Set(times),Set(btotal[!,:t]))...])
        if length(timesMissing) > 0
            for t in timesMissing
                lastTime = maximum(btotal[btotal.t .< t,:t])
                new = btotal[btotal.t .== lastTime,:]
                replace!(new.t,lastTime => t)
                append!(btotal,new)
            end
        end
        bfreq = DataFrame(execute(dbTempSim, "SELECT t, bstrain_id, abundance 
                    FROM babundance 
                    WHERE run_id = $(run_id2) 
                    AND t in ($(join(times,", ")))"))
        rename!(bfreq, :abundance => :bfreq)
        timesMissing = sort([setdiff(Set(times),Set(bfreq[!,:t]))...])
        if length(timesMissing) > 0
            for t in timesMissing
                lastTime = maximum(bfreq[bfreq.t .< t,:t])
                new = bfreq[bfreq.t .== lastTime,:]
                replace!(new.t,lastTime => t)
                append!(bfreq,new)
            end
        end
        bfreq = innerjoin(btotal,bfreq, on=:t)
        bfreq[!,:bfreq] = bfreq[!,:bfreq]./bfreq[!,:btotal]
        select!(bfreq, Not([:btotal]))
        bspacers = DataFrame(execute(dbTempSim, "SELECT bstrain_id, spacer_id 
                    FROM bspacers WHERE run_id = $(run_id2)"))
        for (t,vstrainID,freq) in zip(vfreq[!,:t],vfreq[!,:vstrain_id],vfreq[!,:vfreq])
            if isempty(btotal[btotal.t.==t,:])
                execute(dbOutput, "INSERT INTO local_adaptation VALUES (?,?,?,?)",
                        (run_id2, t, vstrain_id, 0))
                continue
            end
            bspacersNow = innerjoin(bfreq[bfreq.t.==t,:][!,[:t,:bstrain_id]],bspacers, on=:bstrain_id)
            iNet = unique(select(innerjoin(bspacersNow,vspacers[vspacers.vstrain_id .== vstrainID,:], 
                        on=[:spacer_id]), Not([:spacer_id])))
            if isempty(iNet)
                execute(dbOutput, "INSERT INTO local_adaptation VALUES (?,?,?,?)",
                    (run_id2, t, vstrainID, freq))
                continue
            else
                bsus = setdiff(Set(bfreq[bfreq.t.==t,:bstrain_id]),Set(iNet[!,:bstrain_id]))
                if length(bsus) == 0
                    execute(dbOutput, "INSERT INTO local_adaptation VALUES (?,?,?,?)",
                    (run_id2, t, vstrainID, 0))
                    continue
                else
                    fitness = freq*sum(bfreq[(bfreq.t.==t).&(in(bsus).(bfreq.bstrain_id)),:bfreq])
                     execute(dbOutput, "INSERT INTO local_adaptation VALUES (?,?,?,?)",
                    (run_id2, t, vstrainID, fitness))
                end
            end
        end 
    end
end

        # bfreq = DataFrame(execute(dbTempSim, "SELECT t, microbial_abundance 
        #         FROM summary 
        #         WHERE run_id = $(run_id2) AND microbial_abundance != 0
        #         AND t in ($(join(unique(vfreq[!,:t]),", ")))"))
        # if isempty(bfreq)
        #     for (t,vstrainID) in zip(vfreq[!,:t],vfreq[!,:vstrain_id])
        #             execute(dbOutput, "INSERT INTO local_adaptation VALUES (?,?,?,?,?)",
        #             (run_id, run_id2, t, vstrainID, 0))
        #     end
        #     continue
        # end


localAdaptation(run_id,run_ids2,timeInc)
println("Complete!")