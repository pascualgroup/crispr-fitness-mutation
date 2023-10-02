#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
# using Interpolations
# using Dierckx
# using DataInterpolations
using DataFramesMeta

run_id = ARGS[1]
##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbOutputPath = joinpath("epidemic-count_output.sqlite") # cluster
# dbOutputPath = joinpath("/Volumes/Yadgah/walls-shannon_output.sqlite") # local
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/epidemic-count.sqlite") # local

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

function identifyEscape(vstrainID, parentID, tCreate)
    parentHosts = []
    parentSpacers = [sID for (sID,) in
                     execute(dbTempSim,"SELECT spacer_id FROM vpspacers 
                        WHERE run_id = $(run_id) AND vstrain_id = $(parentID)")]
    progHosts = []
    progSpacers = [sID for (sID,) in
                   execute(dbTempSim,"SELECT spacer_id FROM vpspacers 
                        WHERE run_id = $(run_id) AND vstrain_id = $(vstrainID)")]
    for (bstrainID,) in execute(
        dbTempSim,
        "SELECT DISTINCT bstrain_id FROM babundance
        WHERE run_id = $(run_id) AND t = $(tCreate) ORDER BY bstrain_id")
        # println("bstrain: $(bstrain_id)")
        if bstrainID == 1
            continue
        end
        bspacers = [sID for (sID,) in
                    execute(dbTempSim,"SELECT spacer_id FROM bspacers 
            WHERE run_id = $(run_id) AND bstrain_id = $(bstrainID)")]
        if length(intersect(bspacers, parentSpacers)) == 0
            append!(parentHosts, bstrainID)
        end
        if length(intersect(bspacers, progSpacers)) == 0
            append!(progHosts, bstrainID)
        end
    end
    if length(progHosts) > length(parentHosts)
        return true
    else
        return false
    end
end


vseries = DataFrame(execute(dbSim, "SELECT t, viral_abundance FROM summary WHERE run_id = $(run_id) AND viral_abundance != 0"))
vabund = DataFrame(execute(dbSim, "SELECT t, vstrain_id, abundance FROM vabundance WHERE run_id = $(run_id) AND abundance != 0"))
numEscapes = 0
escapes = []

for vstrainID in unique(vabund.vstrain_id)
    # println(vstrainID)
    vstrainAbund = vabund[vabund.vstrain_id.==vstrainID, :][:, :abundance]
    maxAbund = maximum(vstrainAbund)
    maxT = vabund[(vabund.vstrain_id.==vstrainID).&(vabund.abundance.==maxAbund), :][:, :t][1]
    vtotal = vseries[vseries.t.==maxT, :].viral_abundance[1]
    # if maxAbund/vtotal < .005
    #     continue
    # end
    # if maxAbund < 1000
    #     continue
    # end
    (mutant,) = execute(
        dbSim,
        "SELECT t_creation, parent_vstrain_id
        FROM vstrains WHERE run_id = $(run_id) AND vstrain_id = $(vstrainID)"
    )
    tCreate = ceil(mutant.t_creation)
    parentID = mutant.parent_vstrain_id
    escape = identifyEscape(vstrainID, parentID, tCreate)
    # print("vstrainID: $(vstrainID), escape: $(escape)")
    if escape
        # t = vabund[vabund.vstrain_id.==vstrainID, :][peakIdx, :t]
        # println("vstrainID: $(vstrainID), t = $(t)")
        # append!(vlines, t)
        append!(escapes, vstrainID)
        numEscapes += 1
    end
end

println("Complete!")
