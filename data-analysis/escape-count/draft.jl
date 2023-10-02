using SQLite
using DataFrames
using SQLite.DBInterface: execute
using Interpolations
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

vseries = DataFrame(execute(dbSim, "SELECT t, viral_abundance FROM summary WHERE run_id = $(run_id) AND viral_abundance != 0"))
vabund = DataFrame(execute(dbSim, "SELECT t, vstrain_id, abundance FROM vabundance WHERE run_id = $(run_id) AND abundance != 0"))

# scatter(vseries[!, :t], vseries[!, :viral_abundance], label="input data")

plot(fxnDeriv)
plot!(smoothfxn)

smoothfxn = interpolate(vseries.viral_abundance, BSpline(Cubic(Line(OnGrid()))))
totalD = [Interpolations.gradient(smoothfxn, x)[1]
            for x in collect(1:1:length(vabund))]

plot(vseries.viral_abundance)
plot!(vlines, seriestype="vline")


(comboID,) = execute(dbSim, "SELECT combo_id FROM runs WHERE run_id = $(run_id)")
comboID = comboID.combo_id
(beta,) = execute(dbSim, "SELECT viral_burst_size FROM param_combos WHERE combo_id = $(comboID)")
beta = beta.viral_burst_size


function identifyEscape(vstrainID,parentID,tCreate)
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
        if length(intersect(bspacers,parentSpacers)) == 0
            append!(parentHosts,bstrainID)
        end
        if length(intersect(bspacers,progSpacers)) == 0
            append!(progHosts,bstrainID)
        end  
    end
    if length(progHosts) > length(parentHosts)
        return true
    else
        return false
    end
end


vseries = DataFrame(execute(dbTempSim, "SELECT t, viral_abundance FROM summary WHERE run_id = $(run_id) AND viral_abundance != 0"))
vabund = DataFrame(execute(dbTempSim, "SELECT t, vstrain_id, abundance FROM vabundance WHERE run_id = $(run_id) AND abundance != 0"))

# scatter(vseries[!, :t], vseries[!, :viral_abundance], label="input data")

plot(fxnDeriv)
plot!(smoothfxn)
smoothfxn = interpolate(vseries.viral_abundance, BSpline(Cubic(Line(OnGrid()))))
totalD = [Interpolations.gradient(smoothfxn, x)[1]
          for x in collect(1:1:length(vabund))]
plot(vseries.t, vseries.viral_abundance)
plot!(vabund[vabund.vstrain_id.==312, :].t,vabund[vabund.vstrain_id.==312, :].abundance)
plot!(vlines, seriestype="vline")

plot(vabund[vabund.t.==60,:])

(comboID,) = execute(dbSim, "SELECT combo_id FROM runs WHERE run_id = $(run_id)")
comboID = comboID.combo_id
(beta,) = execute(dbSim, "SELECT viral_burst_size FROM param_combos WHERE combo_id = $(comboID)")
beta = beta.viral_burst_size
#
numEpiEscape = 0
numEpi = 0
first = true
vlines = []
vstrains = []
for vstrainID in unique(vabund.vstrain_id)
    # println(vstrainID)
    vstrainAbund = vabund[vabund.vstrain_id.==vstrainID, :][:, :abundance]
    if length(vstrainAbund) < 3
        continue
    end
    smoothfxn = interpolate(vstrainAbund, BSpline(Cubic(Line(OnGrid()))))
    fxnDeriv = [Interpolations.gradient(smoothfxn, x)[1]
                for x in collect(1:1:length(vstrainAbund))]
    boolDeriv = ones(Bool, length(vstrainAbund))
    negIndices = findall(x -> x < 0, fxnDeriv)
    boolDeriv[negIndices] .= false
    epi0 = false
    # diff = vstrainAbund[1:end] .- [vstrainAbund[1], vstrainAbund[1:end-1]...]
    # thresholdIdx = findall(x -> x >= 100, diff)
    # # thresholdIdx = findall(x -> x >= 5, diff)
    # if length(thresholdIdx) == 0
    #     continue
    # else
    #     betaInc = minimum(thresholdIdx)
    # end
    maxAbund = maximum(vstrainAbund)
    maxT = vabund[(vabund.vstrain_id .== vstrainID)\
                .&(vabund.abundance .== maxAbund), :][:, :t][1]
    vtotal = vseries[vseries.t.==maxT, :].viral_abundance
    if maxAbund/vtotal < .01
        continue
    else
    
    firstIdx = findfirst(boolDeriv[betaInc:end])
    if firstIdx !== Nothing
        epi0 = true
    else
        continue
    end
    if epi0
        peakIdx = findfirst(false .== boolDeriv[betaInc:end])
        if peakIdx !== Nothing
            # dt = vabund[vabund.vstrain_id.==vstrainID, :t][1]
            peakIdx = peakIdx + betaInc - 1
            if numEpi > 1
                (mutant,) = execute(
                    dbSim,
                    "SELECT t_creation, parent_vstrain_id
                    FROM vstrains WHERE run_id = $(run_id) AND vstrain_id = $(vstrainID)"
                )
                tCreate = floor(mutant.t_creation)
                parentID = mutant.parent_vstrain_id
                escape = identifyEscape(vstrainID, parentID, tCreate)
                print(escape)
                if escape
                    t = vabund[vabund.vstrain_id.==vstrainID, :][peakIdx, :t]
                    println("vstrainID: $(vstrainID), t = $(t)")
                    append!(vlines, t)
                    append!(vstrains, vstrainID)
                    numEscapes += 1
                end
            end
            numEpi += 1
        end
        continue
    end
end





numEscapes = 0
escapes = []
for vstrainID in unique(vabund.vstrain_id)
    # println(vstrainID)
    vstrainAbund = vabund[vabund.vstrain_id.==vstrainID, :][:, :abundance]
    maxAbund = maximum(vstrainAbund)
    maxT = vabund[(vabund.vstrain_id .== vstrainID).&(vabund.abundance .== maxAbund), :][:, :t][1]
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
    print("vstrainID: $(vstrainID), escape: $(escape)")
    if escape
        # t = vabund[vabund.vstrain_id.==vstrainID, :][peakIdx, :t]
        # println("vstrainID: $(vstrainID), t = $(t)")
        # append!(vlines, t)
        append!(escapes, vstrainID)
        numEscapes += 1
    end
end