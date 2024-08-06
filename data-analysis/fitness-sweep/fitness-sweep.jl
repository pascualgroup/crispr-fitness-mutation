#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]
##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbOutputPath = joinpath("fitness-sweep_output.sqlite") # cluster
rm(dbOutputPath, force=true)
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

# dbRunPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db.sqlite") # cluster
# dbRun = SQLite.DB(dbRunPath)
# (check,) = execute(dbRun, "SELECT combo_id, replicate FROM runs WHERE run_id = $(run_id)")
# dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "runs","c$(check[0])","r$(check[1])","output.sqlite") # cluster
dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/escape-count.sqlite") # local

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
execute(dbOutput, "CREATE TABLE vstrain_escapes (t_creation REAL,  vstrain_id INTEGER, parent_vstrain_id INTEGER, bstrain_id_escaped INTEGER)")
# execute(dbOutput, "CREATE TABLE bstrain_escapes (t_creation REAL,  bstrain_id INTEGER, parent_bstrain_id INTEGER, vstrain_id_escaped INTEGER)")
# execute(dbOutput, "CREATE TABLE extinction_occurrence (microbes INTEGER, viruses INTEGER)")
# execute(dbOutput, "CREATE TABLE simulation_end_time (microbe_end_time REAL, virus_end_time REAL)")
# execute(dbOutput, "CREATE TABLE bstrains (t_creation REAL,  bstrain_id INTEGER, parent_bstrain_id INTEGER, infecting_vstrain_id INTEGER)")
# execute(dbOutput, "CREATE TABLE vstrains (t_creation REAL,  vstrain_id INTEGER, parent_vstrain_id INTEGER, infected_bstrain_id INTEGER)")

function identifyViralEscape(vstrainID, parentID, tCreate)
    parentHosts = []
    parentSpacers = [sID for (sID,) in
                     execute(
        dbSim, "SELECT spacer_id FROM vpspacers WHERE run_id = $(run_id) AND vstrain_id = $(parentID)")]
    progHosts = []
    progSpacers = [sID for (sID,) in
                   execute(
        dbSim,
        "SELECT spacer_id FROM vpspacers WHERE run_id = $(run_id) AND vstrain_id = $(vstrainID)")]
    for (bstrainID,) in execute(
        dbSim,
        "SELECT DISTINCT bstrain_id FROM babundance WHERE run_id = $(run_id) AND t = $(tCreate) ORDER BY bstrain_id")
        # println("bstrain: $(bstrain_id)")
        if bstrainID == 1
            continue
        end
        bspacers = [sID for (sID,) in
                    execute(
            dbSim, "SELECT spacer_id FROM bspacers WHERE run_id = $(run_id) AND bstrain_id = $(bstrainID)")]
        if length(intersect(bspacers, parentSpacers)) == 0
            append!(parentHosts, bstrainID)
        end
        if length(intersect(bspacers, progSpacers)) == 0
            append!(progHosts, bstrainID)
        end
    end
    if length(progHosts) > length(parentHosts)
        return true, setdiff(progHosts, parentHosts)
    else
        return false, nothing
    end
end

function viralEscapes()
    vabund = DataFrame(execute(dbSim, "SELECT DISTINCT vstrain_id FROM vabundance WHERE run_id = $(run_id) AND abundance != 0"))
    numEscapes = 0
    escapes = []
    for vstrainID in unique(vabund.vstrain_id)
        # println(vstrainID)
        # vstrainAbund = vabund[vabund.vstrain_id.==vstrainID, :][:, :abundance]
        # maxAbund = maximum(vstrainAbund)
        # maxT = vabund[(vabund.vstrain_id.==vstrainID).&(vabund.abundance.==maxAbund), :][:, :t][1]
        # vtotal = vseries[vseries.t.==maxT, :].viral_abundance[1]
        # # if maxAbund/vtotal < .005
        # #     continue
        # # end
        # # if maxAbund < 1000
        # #     continue
        # # end
        (mutant,) = execute(
            dbSim,
            "SELECT t_creation, parent_vstrain_id
            FROM vstrains WHERE run_id = $(run_id) AND vstrain_id = $(vstrainID)"
        )
        tCreate = ceil(mutant.t_creation)
        parentID = mutant.parent_vstrain_id
        escape, bstrainsEscaped = identifyViralEscape(vstrainID, parentID, tCreate)
        # print("vstrainID: $(vstrainID), escape: $(escape)")
        if escape
            # t = vabund[vabund.vstrain_id.==vstrainID, :][peakIdx, :t]
            # println("vstrainID: $(vstrainID), t = $(t)")
            # append!(vlines, t)
            append!(escapes, vstrainID)
            numEscapes += 1
            for bstrainID in bstrainsEscaped
                execute(dbOutput, "INSERT INTO vstrain_escapes VALUES (?,?,?,?)",
                    (mutant.t_creation, vstrainID, parentID, bstrainID))
            end
        end
    end
    df = DataFrame(num_escapes=Int64(numEscapes))
    df |> SQLite.load!(dbOutput,
        "viral_escape_count", ifnotexists=true)
end

function identifyHostEscape(bstrainID, parentID, tCreate)
    parentHosts = []
    parentSpacers = [sID for (sID,) in
                     execute(
        dbSim, "SELECT spacer_id FROM bspacers WHERE run_id = $(run_id) AND bstrain_id = $(parentID)")]
    progHosts = []
    progSpacers = [sID for (sID,) in
                   execute(
        dbSim,
        "SELECT spacer_id FROM bspacers WHERE run_id = $(run_id) AND bstrain_id = $(bstrainID)")]
    for (vstrainID,) in execute(
        dbSim,
        "SELECT DISTINCT vstrain_id FROM vabundance WHERE run_id = $(run_id) AND t = $(tCreate) ORDER BY vstrain_id")
        # println("bstrain: $(bstrain_id)")
        vspacers = [sID for (sID,) in
                    execute(
            dbSim, "SELECT spacer_id FROM vpspacers WHERE run_id = $(run_id) AND vstrain_id = $(vstrainID)")]
        if length(intersect(vspacers, parentSpacers)) == 0
            append!(parentHosts, vstrainID)
        end
        if length(intersect(vspacers, progSpacers)) == 0
            append!(progHosts, vstrainID)
        end
    end
    if length(progHosts) < length(parentHosts)
        return true, setdiff(parentHosts, progHosts)
    else
        return false, nothing
    end
end

function hostEscapes()
    babund = DataFrame(execute(dbSim, "SELECT DISTINCT bstrain_id FROM babundance WHERE run_id = $(run_id) AND abundance != 0"))
    numEscapes = 0
    escapes = []
    for bstrainID in unique(babund.bstrain_id)
        (mutant,) = execute(
            dbSim,
            "SELECT t_creation, parent_bstrain_id
            FROM bstrains WHERE run_id = $(run_id) AND bstrain_id = $(bstrainID)"
        )
        tCreate = ceil(mutant.t_creation)
        parentID = mutant.parent_bstrain_id
        escape, vstrainsEscaped = identifyHostEscape(bstrainID, parentID, tCreate)
        if escape
            append!(escapes, bstrainID)
            numEscapes += 1
            for vstrainID in vstrainsEscaped
                execute(dbOutput, "INSERT INTO bstrain_escapes VALUES (?,?,?,?)",
                    (mutant.t_creation, bstrainID, parentID, vstrainID))
            end
        end
    end
    df = DataFrame(num_escapes=Int64(numEscapes))
    df |> SQLite.load!(dbOutput,
        "host_escape_count", ifnotexists=true)
end

# function extinction()
#     microbesDF = DataFrame(execute(
#         dbSim,
#         "SELECT t, microbial_abundance FROM summary WHERE run_id = $(run_id) AND microbial_abundance = 0 ORDER BY t DESC"
#     ))
#     virusesDF = DataFrame(execute(
#         dbSim,
#         "SELECT t, viral_abundance FROM summary WHERE run_id = $(run_id) AND viral_abundance = 0 ORDER by t DESC"
#     ))
#     microbesDF = microbesDF[(microbesDF.t.!=0), :]
#     virusesDF = virusesDF[(virusesDF.t.!=0), :]
#     microbeExt = 0
#     virusExt = 0
#     simEndTime = 0

#     if issubset(0, microbesDF[:, :microbial_abundance])
#         microbeExt = 1
#         # mSimEndTime = maximum(microbesDF[(microbesDF.microbial_abundance.!=0), :].t) + 1
#         println(microbeExt)
#     else
#         microbeExt = 0
#         # mSimEndTime = maximum(microbesDF.t)
#         println(microbeExt)
#     end

#     if issubset(0, virusesDF[:, :viral_abundance])
#         if issubset(0, microbesDF[:, :microbial_abundance])
#             virusExt = 0
#             # vSimEndTime = maximum(virusesDF.t)
#         else
#             virusExt = 1
#             # vSimEndTime = maximum(virusesDF[(virusesDF.viral_abundance.!=0), :].t) + 1
#         end
#         #println(virusExt)
#     else
#         virusExt = 0
#         vSimEndTime = maximum(virusesDF.t)
#         #println(virusExt)
#     end

#     execute(dbOutput, "BEGIN TRANSACTION")
#     execute(dbOutput, "INSERT INTO extinction_occurrence VALUES (?,?)", (microbeExt, virusExt))
#     # execute(dbOutput, "INSERT INTO simulation_end_time VALUES ($(mSimEndTime),$(vSimEndTime))")
#     execute(dbOutput, "COMMIT")

#     tableNames = [table for (table,) in
#                   execute(
#         dbSim,
#         "SELECT name FROM sqlite_master WHERE type='table' AND name in ('bextinctions', 'vextinctions');"
#     )]
#     if length(tableNames) > 0
#         execute(dbOutput, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
#         for table in tableNames
#             if table == "vextinctions"
#                 strainID = "vstrain_id"
#                 execute(dbOutput, "CREATE TABLE $(table) ($(strainID) INTEGER, t_extinction REAL)")
#                 execute(dbOutput, "BEGIN TRANSACTION")
#                 execute(dbOutput, "INSERT INTO $(table)($(strainID), t_extinction) SELECT $(strainID), t_extinction FROM dbSim.$(table) WHERE run_id = $(run_id);")
#                 execute(
#                     dbOutput,
#                     "INSERT INTO vstrains (t_creation,  vstrain_id, parent_vstrain_id, infected_bstrain_id) 
#                     SELECT t_creation,  vstrain_id, parent_vstrain_id, infected_bstrain_id 
#                     FROM dbSim.vstrains WHERE run_id = $(run_id);"
#                 )
#                 execute(dbOutput, "COMMIT")
#             end
#             if table == "bextinctions"
#                 strainID = "bstrain_id"
#                 execute(dbOutput, "CREATE TABLE $(table) ($(strainID) INTEGER, t_extinction REAL)")
#                 execute(dbOutput, "BEGIN TRANSACTION")
#                 execute(dbOutput, "INSERT INTO $(table)($(strainID), t_extinction) SELECT $(strainID), t_extinction FROM dbSim.$(table) WHERE run_id = $(run_id);")
#                 execute(
#                     dbOutput,
#                     "INSERT INTO bstrains (t_creation,  bstrain_id, parent_bstrain_id, infecting_vstrain_id) 
#                     SELECT t_creation,  bstrain_id, parent_bstrain_id, infecting_vstrain_id 
#                     FROM dbSim.bstrains WHERE run_id = $(run_id);"
#                 )
#                 execute(dbOutput, "COMMIT")
#             end
#         end
#     end
# end

# println("Processing extinction occurrences of run $(run_id)")
# extinction()
println("...viral escapes...")
viralEscapes()
# println("...host escapes...")
# hostEscapes()
println("Complete!")
