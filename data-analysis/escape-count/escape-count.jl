#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]
##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbOutputPath = joinpath("escape-count_output.sqlite") # cluster
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/escape-count.sqlite") # local

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
execute(dbOutput, "CREATE TABLE vstrain_escapes (t_creation REAL,  vstrain_id INTEGER, parent_vstrain_id INTEGER, bstrain_id_escaped INTEGER)")
execute(dbOutput, "CREATE TABLE bstrain_escapes (t_creation REAL,  bstrain_id INTEGER, parent_bstrain_id INTEGER, vstrain_id_escaped INTEGER)")

function identifyViralEscape(vstrainID, parentID, tCreate)
    parentHosts = []
    parentSpacers = [sID for (sID,) in
                     execute(
        dbSim, "SELECT spacer_id FROM vpspacers 
         WHERE run_id = $(run_id) AND vstrain_id = $(parentID)")]
    progHosts = []
    progSpacers = [sID for (sID,) in
                   execute(
        dbSim,
        "SELECT spacer_id FROM vpspacers 
        WHERE run_id = $(run_id) AND vstrain_id = $(vstrainID)")]
    for (bstrainID,) in execute(
        dbSim,
        "SELECT DISTINCT bstrain_id FROM babundance
        WHERE run_id = $(run_id) AND t = $(tCreate) ORDER BY bstrain_id")
        # println("bstrain: $(bstrain_id)")
        if bstrainID == 1
            continue
        end
        bspacers = [sID for (sID,) in
                    execute(
            dbSim, "SELECT spacer_id FROM bspacers 
            WHERE run_id = $(run_id) AND bstrain_id = $(bstrainID)")]
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
    # vseries = DataFrame(execute(dbSim, "SELECT t, viral_abundance FROM summary WHERE run_id = $(run_id) AND viral_abundance != 0"))
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
        dbSim, "SELECT spacer_id FROM bspacers 
         WHERE run_id = $(run_id) AND bstrain_id = $(parentID)")]
    progHosts = []
    progSpacers = [sID for (sID,) in
                   execute(
        dbSim,
        "SELECT spacer_id FROM bspacers 
        WHERE run_id = $(run_id) AND bstrain_id = $(bstrainID)")]
    for (vstrainID,) in execute(
        dbSim,
        "SELECT DISTINCT vstrain_id FROM vabundance
        WHERE run_id = $(run_id) AND t = $(tCreate) ORDER BY vstrain_id")
        # println("bstrain: $(bstrain_id)")
        vspacers = [sID for (sID,) in
                    execute(
            dbSim, "SELECT spacer_id FROM vpspacers 
            WHERE run_id = $(run_id) AND vstrain_id = $(vstrainID)")]
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



viralEscapes()
hostEscapes()
println("Complete!")
