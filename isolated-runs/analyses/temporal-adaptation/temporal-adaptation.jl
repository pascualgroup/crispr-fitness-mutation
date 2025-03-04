#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
# using BenchmarkTools


run_id = ARGS[1]
delayInc = parse(Int64, ARGS[2])

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("temporal-adaptation_output.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/9_MOI3", "sweep_db_gathered.sqlite")
# dbTriPath = joinpath("/Volumes/Yadgah/comboID66/tripartite-networksC66.sqlite") # local
# dbShanPath = joinpath("/Volumes/Yadgah/comboID66/shannonC66.sqlite") # local
# dbMatchPath = joinpath("/Volumes/Yadgah/comboID66/matchesC66.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/comboID66/temporal-adaptation.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/comboID1/temporal-adaptation.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/9_MOI3", "temporal-adaptation.sqlite")


#if isfile(dbOutputPath)
#    error("temporal-adaptation_output.sqlite already exists; delete first")
#end # cluster
##
#dbTempSim = SQLite.DB(dbSimPath)
#(cr,) = execute(dbTempSim, "SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
#if length(ARGS) > 1 && ARGS[end] == "combo"
#    if isfile(joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
#        "gathered-analyses", "tripartite-networks", "tripartite-networks.sqlite"))
#        dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
#            "gathered-analyses", "matches", "matches.sqlite")
#        dbTriPath = joinpath(SCRIPT_PATH, "..", "..", "..", "data-analysis",
#            "gathered-analyses", "tripartite-networks", "tripartite-networks.sqlite")
#    elseif isfile(joinpath(SCRIPT_PATH, "..", "..", "isolates",
#        "comboID$(cr[1])", "tripartite-networksC$(cr[1]).sqlite"))
#        dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
#            "comboID$(cr[1])", "matchesC$(cr[1]).sqlite")
#        dbTriPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
#            "comboID$(cr[1])", "tripartite-networksC$(cr[1]).sqlite")
#    else
#        error("Neither /data-analysis/gathered-analyses/tripartite-networks/tripartite-networks.sqlite exists \n
#            nor does /isolated-runs/isolates/comboID$(cr[1])/tripartite-networksC$(cr[1]).sqlite; \n
#            compute one of them first. \n
#            Note that `matches` analysis needs to be run before tripartite network analysis.")
#    end
#end
#if length(ARGS) == 1
#    dbMatchPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
#        "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)", "matches_output.sqlite")
#    dbTriPath = joinpath(SCRIPT_PATH, "..", "..", "isolates",
#        "runID$(run_id)-c$(cr.combo_id)-r$(cr.replicate)", "tripartite-networks_output.sqlite")
#    if !isfile(dbMatchPath)
#        error("matches_output.sqlite does not exist; compute first")
#    end
#
#end
#dbTempMatch = SQLite.DB(dbMatchPath)
dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE temporal_adaptation (time_delay REAL, TA REAL)")
execute(dbOutput, "CREATE TABLE viral_frequency_dependent_fitness (t REAL, vstrain_id INTEGER, fitness REAL)")
execute(dbOutput, "CREATE TABLE bstrain_to_vstrain_0matches (vstrain_id INTEGER, bstrain_id INTEGER)")
maxT = maximum([t for (t,) in execute(dbTempSim, "SELECT t_final FROM param_combos")])

# Create temporary database that is a copy of the main database at the run_id value of the script's argument
tmpPath = joinpath(ENV["SLURM_TMPDIR"], "$(run_id).sqlite")
dbTempSim = SQLite.DB(tmpPath)
# dbTempSim = SQLite.DB("/Volumes/Yadgah/test.sqlite") # local
execute(dbTempSim, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTempSim, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")
execute(dbTempSim, "CREATE TABLE vpspacers (vstrain_id INTEGER, spacer_id INTEGER)")
execute(dbTempSim, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
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

function identify0Matches()
    maxbID = maximum([bID for (bID,) in execute(dbTempSim,
        "SELECT DISTINCT bstrain_id FROM babundance
        ORDER BY bstrain_id")])
    maxsID = maximum([sID for (sID,) in
                      execute(dbTempSim, "SELECT spacer_id FROM vpspacers")])
    maxvID = maximum([bID for (bID,) in execute(dbTempSim,
        "SELECT DISTINCT vstrain_id FROM vabundance
        ORDER BY vstrain_id")])
    vspacers = Array{Array{Bool,1},1}([zeros(Bool, maxsID) for i in 1:maxvID])
    bspacers = Array{Array{Array{Bool,1},1},1}([[zeros(Bool, maxsID) for i in 1:maxbID] for i in 1:maxvID])
    for (bstrain_id,) in execute(
        dbTempSim,
        "SELECT DISTINCT bstrain_id FROM babundance
        ORDER BY bstrain_id")
        # println("bstrain: $(bstrain_id)")
        if bstrain_id == 1
            continue
        end
        for (spacer_id,) in
            execute(dbTempSim, "SELECT spacer_id FROM bspacers WHERE bstrain_id = $(bstrain_id)")
            for (vstrain_id,) in execute(
                dbTempSim,
                "SELECT DISTINCT vstrain_id FROM vabundance
                ORDER BY vstrain_id")
                bspacers[vstrain_id][bstrain_id][spacer_id] = Bool(1)
            end
        end
    end
    for (vstrain_id,) in execute(
        dbTempSim,
        "SELECT DISTINCT vstrain_id FROM vabundance
        ORDER BY vstrain_id")
        # println("vstrain: $(vstrain_id)")
        for (spacer_id,) in
            execute(dbTempSim, "SELECT spacer_id FROM vpspacers WHERE vstrain_id = $(vstrain_id)")
            # println(spacer_id)
            vspacers[vstrain_id][spacer_id] = Bool(1)
        end
    end
    for (vstrain_id,) in execute(
        dbTempSim,
        "SELECT DISTINCT vstrain_id FROM vabundance
        ORDER BY vstrain_id")
        for (bstrain_id,) in execute(
            dbTempSim,
            "SELECT DISTINCT bstrain_id FROM babundance
            ORDER BY bstrain_id")
            # println("vstrain: $(vstrain_id), bstrain: $(bstrain_id)")
            bspacers[vstrain_id][bstrain_id] = bspacers[vstrain_id][bstrain_id] .* vspacers[vstrain_id]
            if sum(bspacers[vstrain_id][bstrain_id]) == 0
                execute(dbOutput, "INSERT INTO bstrain_to_vstrain_0matches VALUES (?,?)",
                    (vstrain_id, bstrain_id))
            end
        end
    end
    execute(dbOutput, "BEGIN TRANSACTION")
    execute(dbOutput, "CREATE INDEX match_index ON bstrain_to_vstrain_0matches (vstrain_id)")
    execute(dbOutput, "COMMIT")
end

function temporalAdaptation(maxT, delayInc)
    maxSimT = maximum([t for (t,) in execute(dbTempSim, "SELECT t FROM vabundance")])
    delays = unique([reverse(collect(0:-delayInc:-maxT))..., collect(0:delayInc:maxT)...])
    for delay in delays
        println("DELAY: $(delay)")
        if maxSimT + 1 - abs(delay) == 0
            continue
        end
        TA = 0
        n = 0
        for (t,) in execute(dbTempSim, "SELECT DISTINCT t FROM vabundance ORDER BY t")
            # println("time: $(t)")
            if t + delay > maxSimT || t + delay < 0
                continue
            end
            btimes = [t for (t,) in execute(dbTempSim, "SELECT t FROM summary")]
            if !in(t+delay,btimes)
                tpast = maximum(btimes[btimes.<t+delay])
                (btotal,) = execute(
                dbTempSim,
                "SELECT microbial_abundance FROM summary
                WHERE t = $(tpast)"
                )
            else
                (btotal,) = execute(
                dbTempSim,
                "SELECT microbial_abundance FROM summary
                WHERE t = $(t + delay)"
                )
            end
            btotal = btotal.microbial_abundance
            (vtotal,) = execute(
                dbTempSim,
                "SELECT viral_abundance FROM summary
                WHERE t = $(t)"
            )
            vtotal = vtotal.viral_abundance
            if vtotal == 0 || btotal == 0
                continue
            end 
            vstrains = DataFrame(execute(
                dbTempSim,
                "SELECT vstrain_id, abundance FROM vabundance
                WHERE t = $(t)"
            ))
            rename!(vstrains, :abundance => :vfreq)
            vstrains = vstrains[vstrains.vfreq.!=0, :]
            if !in(t+delay,btimes)
                tpast = maximum(btimes[btimes.<t+delay])
                bstrains = DataFrame(execute(
                    dbTempSim,
                    "SELECT bstrain_id, abundance FROM babundance
                    WHERE t = $(tpast)"
                ))
                rename!(bstrains, :abundance => :bfreq)
                bstrains = bstrains[bstrains.bfreq.!=0, :]
                (btotal,) = execute(
                    dbTempSim,
                    "SELECT microbial_abundance FROM summary
                    WHERE t = $(tpast)"
                )
            else
                bstrains = DataFrame(execute(
                    dbTempSim,
                    "SELECT bstrain_id, abundance FROM babundance
                    WHERE t = $(t + delay)"
                ))
                rename!(bstrains, :abundance => :bfreq)
                bstrains = bstrains[bstrains.bfreq.!=0, :]
                (btotal,) = execute(
                    dbTempSim,
                    "SELECT microbial_abundance FROM summary
                    WHERE t = $(t + delay)"
                )
            end
            btotal = btotal.microbial_abundance
            matches = DataFrame(execute(
                dbOutput,
                "SELECT bstrain_id, vstrain_id FROM bstrain_to_vstrain_0matches
                WHERE vstrain_id in ($(join(vstrains[!,:vstrain_id],", ")))
                AND bstrain_id in ($(join(bstrains[!,:bstrain_id],", ")))"
            ))
            matches = innerjoin(matches, vstrains, on=:vstrain_id)
            matches = innerjoin(matches, bstrains, on=:bstrain_id)
            if isempty(matches)
                continue
            end
            matches[!, :vfreq] = matches[!, :vfreq] / vtotal
            matches[!, :bfreq] = matches[!, :bfreq] / btotal
            matches[!, :fitness] = matches[!, :bfreq] .* matches[!, :vfreq]
            if delay == 0
                matches = groupby(matches, :vstrain_id)
                matches = combine(matches, [:fitness] .=> sum; renamecols=false)
                matches[!, :t] = fill(Float64(t),length(matches[!,:vstrain_id]))
                matches |> SQLite.load!(dbOutput, "viral_frequency_dependent_fitness", ifnotexists=true)
            end
            TA += sum(matches[!, :fitness])
            n += 1
        end
        if n > 0
            execute(dbOutput, "INSERT INTO temporal_adaptation VALUES (?,?)",
                (delay, TA / (maxSimT + 1 - abs(delay))))
        end
    end
end

identify0Matches()
temporalAdaptation(maxT, delayInc)
# execute(dbOutput, "DROP TABLE bstrain_to_vstrain_0matches")
rm(tmpPath, force=true)
println("Complete!")