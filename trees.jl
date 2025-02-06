#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using Combinatorics


combo_id = parse(Int64,ARGS[1])
dbSimDir = ARGS[2]
bthreshold = parse(Float64,ARGS[3])/100
vthreshold = parse(Float64,ARGS[4])/100

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,dbSimDir,"sweep_db_gathered.sqlite") # cluster
dbSim = SQLite.DB(dbSimPath)

run_id = rand([runID for (runID,) in execute(dbSim,"SELECT run_id FROM runs WHERE combo_id = $(combo_id)")])
run_id = 4314#4526
println("runID is $(run_id)")

dbOutputPath = joinpath(SCRIPT_PATH, "$(dbSimDir)", "trees_output_cID$(combo_id)-runID$(run_id)-bthresh$(ARGS[3])-vthresh$(ARGS[4]).sqlite") # cluster
if isfile(dbOutputPath)
    rm(dbOutputPath,force=true)
end
dbOutput = SQLite.DB(dbOutputPath)
execute(dbOutput, "CREATE TABLE tree_babundance (t REAL, tree_bstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE tree_vabundance (t REAL, tree_vstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE bstrain_creation_extinction (bstrain_id INTEGER, t_creation REAL, t_extinction REAL)")
# execute(dbOutput, "CREATE TABLE tree_bstrain_creation_extinction
# (tree_bstrain_id INTEGER, t_creation REAL, t_extinction REAL,
# tree_parent_bstrain_id INTEGER, tree_infecting_vstrain_id)")
execute(dbOutput, "CREATE TABLE tree_bstrain_creation_extinction
(tree_bstrain_id INTEGER, t_creation REAL, t_extinction REAL,
tree_parent_bstrain_id INTEGER)")
execute(dbOutput, "CREATE TABLE vstrain_creation_extinction (vstrain_id INTEGER, t_creation REAL, t_extinction REAL)")
execute(dbOutput, "CREATE TABLE tree_vstrain_creation_extinction
(tree_vstrain_id INTEGER, t_creation REAL, t_extinction REAL,
tree_parent_vstrain_id INTEGER, tree_infected_bstrain_id INTEGER)")
# execute(dbOutput, "CREATE TABLE MRCA_bstrains
#                     (t REAL, bstrain_id INTEGER, lineage INTEGER, MRCA_bstrain_id INTEGER,
#                         t_creation REAL)")
# execute(dbOutput, "CREATE TABLE MRCA_vstrains
#                     (t REAL, vstrain_id INTEGER, lineage INTEGER, MRCA_vstrain_id INTEGER,
#                         t_creation REAL)")
# execute(dbOutput, "CREATE TABLE pairwise_MRCA_bstrains
#                     (t REAL, bstrain_id_sample_1 INTEGER, bstrain_id_sample_2 INTEGER,
#                     lineage INTEGER, MRCA_bstrain_id INTEGER, t_creation REAL)")
# execute(dbOutput, "CREATE TABLE pairwise_MRCA_vstrains
#                     (t REAL, vstrain_id_sample_1 INTEGER, vstrain_id_sample_2 INTEGER,
#                     lineage INTEGER, MRCA_vstrain_id INTEGER, t_creation REAL)")



# Create temporary database that is a copy of the main database at the run_id value of the script's argument
dbTemp = SQLite.DB()
# rm("/Volumes/Yadgah/test.sqlite",force=true)
# dbTemp = SQLite.DB("/Volumes/Yadgah/test.sqlite") # local
# dbTemp = SQLite.DB(dbSimPath) # local
# execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
# execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")
execute(dbTemp, "CREATE TABLE bstrains (t_creation REAL, bstrain_id INTEGER, parent_bstrain_id INTEGER, infecting_vstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE vstrains (t_creation REAL, vstrain_id INTEGER, parent_vstrain_id INTEGER, infected_bstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE bgrowthrates (bstrain_id INTEGER, growth_rate REAL)")
execute(dbTemp, "CREATE TABLE bextinctions (bstrain_id INTEGER, t_extinction REAL)")
execute(dbTemp, "CREATE TABLE vextinctions (vstrain_id INTEGER, t_extinction REAL)")
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE bspacers (bstrain_id INTEGER, spacer_id INTEGER)")
execute(dbTemp, "CREATE TABLE vpspacers (vstrain_id INTEGER, spacer_id INTEGER)")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
# execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
# execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO summary (t, microbial_abundance, viral_abundance)
SELECT t, microbial_abundance, viral_abundance
FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bstrains (t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id)
SELECT t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id
FROM dbSim.bstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vstrains (t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id)
SELECT t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id
FROM dbSim.vstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bgrowthrates (bstrain_id, growth_rate)
SELECT bstrain_id, growth_rate
FROM dbSim.bgrowthrates WHERE run_id = $(run_id);")
# execute(dbTemp,"INSERT INTO bextinctions (bstrain_id, t_extinction)
# SELECT bstrain_id, t_extinction
# FROM dbSim.bextinctions WHERE run_id = $(run_id);")
# execute(dbTemp,"INSERT INTO vextinctions (vstrain_id, t_extinction)
# SELECT vstrain_id, t_extinction
# FROM dbSim.vextinctions WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance)
SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance)
SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bspacers (bstrain_id, spacer_id)
SELECT bstrain_id, spacer_id FROM dbSim.bspacers WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vpspacers (vstrain_id, spacer_id)
SELECT vstrain_id, spacer_id FROM dbSim.vpspacers WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")
execute(dbTemp, "DETACH DATABASE dbSim")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX summary_index ON summary (t)")
execute(dbTemp, "CREATE INDEX babundance_index ON babundance (t,bstrain_id)")
execute(dbTemp, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTemp, "CREATE INDEX bstrains_index ON bstrains (t_creation, bstrain_id,parent_bstrain_id)")
execute(dbTemp, "CREATE INDEX vstrains_index ON vstrains (t_creation, vstrain_id,parent_vstrain_id)")
execute(dbTemp, "CREATE INDEX bgrowth_index ON bgrowthrates (bstrain_id)")
execute(dbTemp, "COMMIT")


function keepWhichStrains(bthreshold::Float64,vthreshold::Float64)
    keepVstrains = innerjoin(DataFrame(execute(dbTemp, "SELECT t, vstrain_id, abundance FROM vabundance")),
    DataFrame(execute(dbTemp, "SELECT t, viral_abundance FROM summary")), on=:t)
    select!(keepVstrains, Not([:t]))
    rename!(keepVstrains, :abundance => :freq)
    keepVstrains[!, :freq] = keepVstrains[:, :freq] ./ keepVstrains[:, :viral_abundance]
    select!(keepVstrains, Not([:viral_abundance]))
    keepVstrains = combine(groupby(keepVstrains, [:vstrain_id]), :freq .=> maximum, renamecols=false)
    keepVstrains = keepVstrains[keepVstrains[!, :freq].>=vthreshold, :]
    keepVstrains = keepVstrains[:, :vstrain_id]
    parentTrack = copy(keepVstrains)
    while length(parentTrack) != 0
        parentTrack = [strainID for (strainID,) in 
                        execute(dbTemp,"SELECT parent_vstrain_id 
                        FROM vstrains 
                        WHERE vstrain_id in ($(join(parentTrack,", ")))")]
        parentTrack = parentTrack[parentTrack.!=0]
        append!(keepVstrains,parentTrack)
        unique!(keepVstrains)  
    end
    ##
    infectedBstrains = [strainID for (strainID,) in execute(dbTemp, "SELECT infected_bstrain_id 
                                    FROM vstrains WHERE vstrain_id in ($(join(keepVstrains,", ")))")]
    infectedBstrains = unique(infectedBstrains[infectedBstrains.!=0])
    ##
    freq = innerjoin(DataFrame(execute(dbTemp, "SELECT t, bstrain_id, abundance FROM babundance")),
    DataFrame(execute(dbTemp, "SELECT t, microbial_abundance FROM summary")), on=:t)
    select!(freq, Not([:t]))
    rename!(freq,:abundance=>:freq)
    freq[!,:freq] = freq[:,:freq]./freq[:,:microbial_abundance]
    select!(freq, Not([:microbial_abundance]))
    freq = combine(groupby(freq,[:bstrain_id]),:freq.=>maximum, renamecols=false)
    freq = freq[freq[!,:freq].>=bthreshold,:]
    keepStrains = freq[:,:bstrain_id]
    ##
    append!(keepStrains, infectedBstrains)
    unique!(keepStrains)
    ##
    parentTrack = copy(keepStrains)
    while length(parentTrack) != 0
        parentTrack = [strainID for (strainID,) in 
                        execute(dbTemp,"SELECT parent_bstrain_id 
                        FROM bstrains 
                        WHERE bstrain_id in ($(join(parentTrack,", ")))")]
        parentTrack = parentTrack[parentTrack.!=0]
        append!(keepStrains,parentTrack)
        unique!(keepStrains)  
    end
    return keepStrains, keepVstrains
end

function growthOrder(initialStrains::Array{Int64,1}, strainType::String, keepStrains::Union{Array{Int64,1},Nothing}, keepVstrains::Union{Array{Int64,1},Nothing})
    if strainType == "virus"
        initialBStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT bstrain_id
                            FROM bstrains 
                            WHERE parent_bstrain_id = 0 
                            ORDER BY bstrain_id")
                        ]
        initialBStrainsGrowthOrdered =[strain_id for (strain_id,) in execute(dbTemp,"SELECT bstrain_id
                        FROM bgrowthrates 
                        WHERE bstrain_id in ($(join(initialBStrains,", "))) 
                        ORDER BY growth_rate DESC")
                    ]
        spacers = [[spacerID for (spacerID,) in execute(dbTemp,"SELECT spacer_id
                            FROM bspacers 
                            WHERE bstrain_id = $(bstrainID)")] 
                                for bstrainID in initialBStrainsGrowthOrdered
                        ]
        # noSpacersIdx = findall(x -> x == [], spacers)
        immune =  [[vstrainID for (vstrainID,) in execute(dbTemp,"SELECT vstrain_id
                            FROM vpspacers 
                            WHERE vstrain_id in ($(join(initialStrains,", ")))
                            AND spacer_id in ($(join(spacerSet,", ")))")] 
                                for spacerSet in spacers
                        ]
                              
        initialStrains = unique(reduce(vcat,[setdiff(initialStrains,immuneStrains) for immuneStrains in immune]))
        return initialStrains
    end
    if strainType == "microbe"
        initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT bstrain_id
                        FROM bgrowthrates 
                        WHERE bstrain_id in ($(join(initialStrains,", "))) 
                        ORDER BY growth_rate DESC")
                        ]
        println("initialStrainGrowthOrder: $(initialStrains)")
        return initialStrains
    end
end

function orderTreeStrains(strainType::String, keepStrains::Union{Array{Int64,1}, Nothing}, keepVstrains::Array{Int64,1})
    if strainType == "virus"
        s = "v"
        println("generating viral tree...")
        initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
                        FROM $(s)strains 
                        WHERE parent_$(s)strain_id = 0 
                        ORDER BY $(s)strain_id")
                    ]
        initialStrains = growthOrder(initialStrains, strainType, keepStrains, keepVstrains)
        strainTreeOrder = DataFrame(vstrain_id = Int64[0,1], tree_vstrain_id = Int64[0,1])
        lineages = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
                                    WHERE $(s)strain_id in ($(join(keepVstrains,", ")))
                                    ORDER BY t_creation")
            lineages[strain_id] = searchLineage(strain_id,strainType,keepVstrains)
        end
        initialtreestrainIDs = Int64[0, 1]
        # for initialStrain in initialStrains
        for strainIDX in collect(2:1:length(initialStrains))
            # if initialStrain == 1
            #     continue
            # end
            strainBefore = initialStrains[strainIDX-1]
            newTreeID = initialtreestrainIDs[strainIDX] + length(lineages[strainBefore])
            push!(initialtreestrainIDs,newTreeID)
        end
        strainTreeOrder = DataFrame(vstrain_id = Int64[0,initialStrains...],
                                        tree_vstrain_id = initialtreestrainIDs)
        for initialStrain in initialStrains
            println("...searching descendents of strain $(initialStrain)...")
            MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
                            FROM $(s)strains 
                            WHERE parent_$(s)strain_id = $(initialStrain)
                            AND $(s)strain_id in ($(join(keepVstrains,", "))) 
                            ORDER BY $(s)strain_id")
                        ]
            while length(MRCA_ids) > 0
                println("...finding positions for strains: $(MRCA_ids)...")
                nextGeneration!(strainTreeOrder,MRCA_ids,lineages,strainType,keepVstrains)
                MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
                            FROM $(s)strains 
                            WHERE parent_$(s)strain_id in ($(join(MRCA_ids,", ")))
                            AND $(s)strain_id in ($(join(keepVstrains,", "))) 
                            ORDER BY $(s)strain_id")
                            ]
                if length(MRCA_ids) == 0
                    println("...all strain positions found...")
                end
            end
        end
        strainTreeOrder |> SQLite.load!(dbOutput,"tree_vstrain_order",ifnotexists=true)
        newAbundanceTimeSeries(strainType, keepVstrains)
        # MRCAtimeSeries(lineages,strainType)
    end
    if strainType == "microbe"
        s = "b"
        ##
        println("generating microbial tree...")
        initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
                            FROM $(s)strains 
                            WHERE parent_$(s)strain_id = 0 
                            ORDER BY $(s)strain_id")
                        ]
        initialStrains = growthOrder(initialStrains, strainType, keepStrains, nothing)
        println("initialStrainGrowthOrder: $(initialStrains)")
        strainTreeOrder = DataFrame(bstrain_id = Int64[0,1], tree_bstrain_id = Int64[0,1])
        lineages = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains 
                                    WHERE $(s)strain_id in ($(join(keepStrains,", ")))
                                    ORDER BY t_creation")
            lineages[strain_id] = searchLineage(strain_id,strainType,keepStrains)
        end
        initialtreestrainIDs = Int64[0, 1]
        # for initialStrain in initialStrains
        for strainIDX in collect(2:1:length(initialStrains))
            # if initialStrain == 1
            #     continue
            # end
            strainBefore = initialStrains[strainIDX-1]
            newTreeID = initialtreestrainIDs[strainIDX] + length(lineages[strainBefore])
            push!(initialtreestrainIDs,newTreeID)
        end
        strainTreeOrder = DataFrame(bstrain_id = Int64[0,initialStrains...],
                                        tree_bstrain_id = initialtreestrainIDs)
        for initialStrain in initialStrains
            println("...searching descendents of strain $(initialStrain)...")
            MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
                            FROM $(s)strains 
                            WHERE parent_$(s)strain_id = $(initialStrain)
                            AND $(s)strain_id in ($(join(keepStrains,", "))) 
                            ORDER BY $(s)strain_id")
                        ]
            while length(MRCA_ids) > 0
                println("...finding positions for strains: $(MRCA_ids)...")
                nextGeneration!(strainTreeOrder,MRCA_ids,lineages,strainType,keepStrains)
                MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
                            FROM $(s)strains 
                            WHERE parent_$(s)strain_id in ($(join(MRCA_ids,", ")))
                            AND $(s)strain_id in ($(join(keepStrains,", "))) 
                            ORDER BY $(s)strain_id")
                            ]
                if length(MRCA_ids) == 0
                    println("...all strain positions found...")
                end
            end
        end
        strainTreeOrder |> SQLite.load!(dbOutput,"tree_bstrain_order",ifnotexists=true)
        newAbundanceTimeSeries(strainType,keepStrains)
        # MRCAtimeSeries(lineages,strainType)
    end
end


function nextGeneration!(strainTreeOrder::DataFrame,
    MRCA_ids::Array{Int64,1},lineages::Dict,strainType::String,
    keepStrains::Array{Int64,1})
    # MRCA_ids = sort([Int64(key) for key in keys(lineages)])
    if strainType == "virus"
        s = "v"
        parentStrainIDs = Int64[]
        for MRCA_id in MRCA_ids
            (strain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(MRCA_id) 
            AND $(s)strain_id in ($(join(keepStrains,", "))) ORDER BY parent_$(s)strain_id")
            union!(parentStrainIDs,strain.parent_vstrain_id)
        end
        sort!(parentStrainIDs)
        # println("Parents of MRCA are: $(parentStrainIDs)")
        treePositions!(strainTreeOrder, parentStrainIDs, lineages, strainType, keepStrains)
    end
    if strainType == "microbe"
        s = "b"
        parentStrainIDs = Int64[]
        for MRCA_id in MRCA_ids
            (strain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(MRCA_id) 
            AND $(s)strain_id in ($(join(keepStrains,", "))) ORDER BY parent_$(s)strain_id")
            union!(parentStrainIDs,strain.parent_bstrain_id)
        end
        sort!(parentStrainIDs)
        # println("Parents of MRCA are: $(parentStrainIDs)")
        treePositions!(strainTreeOrder, parentStrainIDs, lineages, strainType, keepStrains)
    end
end

function treePositions!(strainTreeOrder::DataFrame,
    parentStrainIDs::Array{Int64,1}, lineages::Dict, strainType::String, 
    keepStrains::Array{Int64,1})
    if strainType == "virus"
        s = "v"
        for parentID in parentStrainIDs
            # println("Locating descendants of strain $(parentID)")
            subMRCA_ids = [strain for (strain,) in execute(dbTemp, "SELECT $(s)strain_id FROM $(s)strains
            WHERE parent_$(s)strain_id = $(parentID) 
            AND $(s)strain_id in ($(join(keepStrains,", "))) 
            ORDER BY $(s)strain_id")]
            firstTreeStrainID = strainTreeOrder[strainTreeOrder.vstrain_id .== parentID,:].tree_vstrain_id[1]
            for subMRCA in subMRCA_ids
                if subMRCA == last(subMRCA_ids)
                    push!(strainTreeOrder,[subMRCA, firstTreeStrainID + 1])
                    continue
                end
                treeStrainID = 1 + firstTreeStrainID + sum([length(lineages[nestMRCA_id])
                for nestMRCA_id in subMRCA_ids if nestMRCA_id > subMRCA])
                push!(strainTreeOrder,[subMRCA, treeStrainID])
            end
        end
    end
    if strainType == "microbe"
        s = "b"
        for parentID in parentStrainIDs
            # println("Locating descendants of strain $(parentID)")
            subMRCA_ids = [strain for (strain,) in execute(dbTemp, "SELECT $(s)strain_id FROM $(s)strains
            WHERE parent_$(s)strain_id = $(parentID) 
            AND $(s)strain_id in ($(join(keepStrains,", "))) 
            ORDER BY $(s)strain_id")]
            firstTreeStrainID = strainTreeOrder[strainTreeOrder.bstrain_id .== parentID,:].tree_bstrain_id[1]
            for subMRCA in subMRCA_ids
                if subMRCA == last(subMRCA_ids)
                    push!(strainTreeOrder,[subMRCA, firstTreeStrainID + 1])
                    continue
                end
                treeStrainID = 1 + firstTreeStrainID + sum([length(lineages[nestMRCA_id])
                for nestMRCA_id in subMRCA_ids if nestMRCA_id > subMRCA])
                push!(strainTreeOrder,[subMRCA, treeStrainID])
            end
        end
    end
end

function newAbundanceTimeSeries(strainType::String, keepStrains::Array{Int64,1})
    if strainType == "virus"
        s = "v"
        for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")
            # println("Changing strain IDs at time = $(time)")
            for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)abundance
                WHERE t = $(time) 
                AND $(s)strain_id in ($(join(keepStrains,", "))) 
                ORDER BY $(s)strain_id")
                (new,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
                (abundance,) = execute(dbTemp, "SELECT abundance
                FROM $(s)abundance WHERE $(s)strain_id = $(strain) AND t = $(time)")
                execute(dbOutput, "INSERT INTO tree_$(s)abundance VALUES (?,?,?)",
                (time,new.tree_vstrain_id,abundance.abundance))
            end
        end
    end
    if strainType == "microbe"
        s = "b"
        for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")
            # println("Changing strain IDs at time = $(time)")
            for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)abundance
                WHERE t = $(time) 
                AND $(s)strain_id in ($(join(keepStrains,", "))) 
                ORDER BY $(s)strain_id")
                (new,) = execute(dbOutput, "SELECT tree_$(s)strain_id
                FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
                (abundance,) = execute(dbTemp, "SELECT abundance
                FROM $(s)abundance WHERE $(s)strain_id = $(strain) AND t = $(time)")
                execute(dbOutput, "INSERT INTO tree_$(s)abundance VALUES (?,?,?)",
                (time,new.tree_bstrain_id,abundance.abundance))
            end
        end
    end
end

# function findExtinctionTimes(strainType::String)
#     if strainType == "virus"
#         s = "v"
#         for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)extinctions
#                                     ORDER BY $(s)strain_id")
#             (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
#             WHERE $(s)strain_id = $(strain)")
#             (tExt,) = execute(dbTemp, "SELECT t_extinction FROM $(s)extinctions
#             WHERE $(s)strain_id = $(strain)")
#             execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
#             VALUES (?,?,?)",(strain,t.t_creation,tExt.t_extinction))
#             (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
#             WHERE $(s)strain_id = $(strain)")
#             (ibstrain,) = execute(dbTemp, "SELECT infected_bstrain_id FROM $(s)strains
#             WHERE $(s)strain_id = $(strain)")
#             (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
#             FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_vstrain_id)")
#             (ibstrain,) = execute(dbOutput, "SELECT tree_bstrain_id
#             FROM tree_bstrain_order WHERE bstrain_id = $(ibstrain.infected_bstrain_id)")
#             (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
#             FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
#             execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
#             VALUES (?,?,?,?,?)",(tree.tree_vstrain_id,t.t_creation,tExt.t_extinction,
#             pstrain.tree_vstrain_id,ibstrain.tree_bstrain_id))
#         end
#     end
#     if strainType == "microbe"
#         s = "b"
#         println("Searching for strain extinctions")
#         for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)extinctions 
#                                     ORDER BY $(s)strain_id")
#             (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
#             WHERE $(s)strain_id = $(strain)")
#             (tExt,) = execute(dbTemp, "SELECT t_extinction FROM $(s)extinctions
#             WHERE $(s)strain_id = $(strain)")
#             execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
#             VALUES (?,?,?)",(strain,t.t_creation,tExt.t_extinction))
#             (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
#             WHERE $(s)strain_id = $(strain)")
#             (ivstrain,) = execute(dbTemp, "SELECT infecting_vstrain_id FROM $(s)strains
#             WHERE $(s)strain_id = $(strain)")
#             (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
#             FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_bstrain_id)")
#             (ivstrain,) = execute(dbOutput, "SELECT tree_vstrain_id
#             FROM tree_vstrain_order WHERE vstrain_id = $(ivstrain.infecting_vstrain_id)")
#             (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
#             FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
#             execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
#             VALUES (?,?,?,?,?)",(tree.tree_bstrain_id,t.t_creation,tExt.t_extinction,
#             pstrain.tree_bstrain_id,ivstrain.tree_vstrain_id))
#         end
#     end
# end

function findExtinctionTimes(strainType::String,keepStrains::Array{Int64,1})
    if strainType == "virus"
        s = "v"
        for (strain,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)abundance
                                    WHERE $(s)strain_id in ($(join(keepStrains,", "))) 
                                    ORDER BY $(s)strain_id")
            (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            tExt = maximum([t for (t,) in execute(dbTemp, "SELECT t FROM $(s)abundance
            WHERE $(s)strain_id = $(strain)")])
            execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
            VALUES (?,?,?)",(strain,t.t_creation,tExt))
            (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            (ibstrain,) = execute(dbTemp, "SELECT infected_bstrain_id FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
            FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_vstrain_id)")
            (ibstrain,) = execute(dbOutput, "SELECT tree_bstrain_id
            FROM tree_bstrain_order WHERE bstrain_id = $(ibstrain.infected_bstrain_id)")
            (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
            FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
            execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
            VALUES (?,?,?,?,?)",(tree.tree_vstrain_id,t.t_creation,tExt,
            pstrain.tree_vstrain_id,ibstrain.tree_bstrain_id))
        end
    end
    if strainType == "microbe"
        s = "b"
        println("Searching for strain extinctions")
        for (strain,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)abundance 
                                    WHERE $(s)strain_id in ($(join(keepStrains,", "))) 
                                    ORDER BY $(s)strain_id")
            # println(strain)
            (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            tExt = maximum([t for (t,) in execute(dbTemp, "SELECT t FROM $(s)abundance
            WHERE $(s)strain_id = $(strain)")])
            execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
            VALUES (?,?,?)",(strain,t.t_creation,tExt))
            (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            # (ivstrain,) = execute(dbTemp, "SELECT infecting_vstrain_id FROM $(s)strains
            # WHERE $(s)strain_id = $(strain)")
            (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
            FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_bstrain_id)")
            # (ivstrain,) = execute(dbOutput, "SELECT tree_vstrain_id
            # FROM tree_vstrain_order WHERE vstrain_id = $(ivstrain.infecting_vstrain_id)")
            (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
            FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
            # execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
            # VALUES (?,?,?,?,?)",(tree.tree_bstrain_id,t.t_creation,tExt,
            # pstrain.tree_bstrain_id,ivstrain.tree_vstrain_id))
            execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
            VALUES (?,?,?,?)",(tree.tree_bstrain_id,t.t_creation,tExt,
            pstrain.tree_bstrain_id))
        end
    end
end

function searchLineage(MRCA_id::Int64,strainType::String,keepStrains)
    if strainType == "virus"
        type = "v"
    end
    if strainType == "microbe"
        type = "b"
    end
    numDescendents = 1
    lineage = Int64[MRCA_id]
    ancestors = Int64[MRCA_id]
    descendents = Int64[]
    while numDescendents > 0
        for id in ancestors
            desc_ids = [id for (id,) in
            execute(dbTemp,"SELECT $(type)strain_id FROM $(type)strains
            WHERE parent_$(type)strain_id in ($(id)) 
            AND $(type)strain_id in ($(join(keepStrains,", "))) ORDER BY $(type)strain_id")]
            append!(descendents,desc_ids)
            lineage = append!(lineage,desc_ids)
            unique!(lineage)
        end
        numDescendents = length(descendents)
        ancestors = descendents[:]
        descendents = Int64[]
    end
    return sort(lineage)
end

# function searchAncestry(daughter_id::Int64,strainType::String)
#     if strainType == "virus"
#         type = "v"
#     end
#     if strainType == "microbe"
#         type = "b"
#     end
#     ancestors = Int64[daughter_id]
#     ancestor_id = daughter_id
#     while ancestor_id > 0
#         ancestor_id = [id for (id,) in
#         execute(dbTemp,"SELECT parent_$(type)strain_id FROM $(type)strains
#         WHERE $(type)strain_id in ($(ancestor_id)) 
#         AND $(type)strain_id in ($(join(keepStrains,", "))) ORDER BY $(type)strain_id")][1]
#         push!(ancestors,ancestor_id)
#     end
#     return sort(ancestors)
# end

# function MRCAtimeSeries(lineages::Dict{Int64,Array{Int64,1}},strainType::String)
#     if strainType == "virus"
#         s = "v"
#         ancestry = Dict{Int64,Array{Int64,1}}()
#         for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
#             ORDER BY t_creation")
#             anc = searchAncestry(strain_id,strainType)
#             ancestry[strain_id] = anc[anc.!=strain_id]
#         end
#         initialStrains = [strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
#         FROM $(s)strains WHERE parent_$(s)strain_id = 0 ORDER BY $(s)strain_id")]
#         for (time,) in execute(dbTemp,"SELECT DISTINCT t FROM $(s)abundance ORDER BY t")
#             strains = [strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
#             FROM $(s)abundance WHERE t = $(time) ORDER BY $(s)strain_id")]
#             if length(strains) == 0
#                 continue
#             end
#             ancestors = map(x->ancestry[x],strains)
#             for lineage in initialStrains
#                 strain_idx = findall(x-> lineage in x, ancestors)
#                 daughters = strains[strain_idx]
#                 MRCA = ancestors[strain_idx]
#                 if length(MRCA) > 1
#                     MRCA = intersect(MRCA...)
#                     MRCA = MRCA[MRCA .!= [0]]
#                     MRCA = maximum(MRCA)
#                 else
#                     MRCA = lineage
#                 end
#                 (t_creation,) = execute(dbTemp,"SELECT t_creation FROM $(s)strains
#                                         WHERE $(s)strain_id = $(MRCA)")
#                 execute(dbOutput, "BEGIN TRANSACTION")
#                 for daughter in daughters
#                     execute(dbOutput, "INSERT INTO MRCA_$(s)strains VALUES (?,?,?,?,?)",
#                                             (time,daughter,lineage,MRCA,t_creation.t_creation))
#                 end
#                 strainpairs = collect(combinations(strain_idx,2))
#                 for pair in strainpairs
#                     strainID1 = pair[1]
#                     strainID2 = pair[2]
#                     daughters = strains[pair]
#                     MRCA = ancestors[[strainID1,strainID2]]
#                     if length(MRCA) > 1
#                         MRCA = intersect(MRCA...)
#                         MRCA = MRCA[MRCA .!= [0]]
#                         MRCA = maximum(MRCA)
#                     else
#                         MRCA = lineage
#                     end
#                     (t_creation,) = execute(dbTemp,"SELECT t_creation FROM $(s)strains
#                                             WHERE $(s)strain_id = $(MRCA)")
#                     execute(dbOutput, "INSERT INTO pairwise_MRCA_$(s)strains VALUES (?,?,?,?,?,?)",
#                             (time,daughters[1],daughters[2],
#                             lineage,MRCA,t_creation.t_creation))
#                 end
#                 execute(dbOutput, "COMMIT")
#             end
#         end
#     end
#     if strainType == "microbe"
#         s = "b"
#         ancestry = Dict{Int64,Array{Int64,1}}()
#         for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
#             ORDER BY t_creation")
#             anc = searchAncestry(strain_id,strainType)
#             ancestry[strain_id] = anc[anc.!=strain_id]
#         end
#         initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
#         FROM $(s)strains WHERE parent_$(s)strain_id = 0 ORDER BY $(s)strain_id")]
#         for (time,) in execute(dbTemp,"SELECT DISTINCT t FROM $(s)abundance ORDER BY t")
#             strains = [strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
#             FROM $(s)abundance WHERE t = $(time) ORDER BY $(s)strain_id")]
#             if length(strains) == 0
#                 continue
#             end
#             ancestors = map(x->ancestry[x],strains)
#             for lineage in initialStrains
#                 strain_idx = findall(x-> lineage in x, ancestors)
#                 daughters = strains[strain_idx]
#                 MRCA = ancestors[strain_idx]
#                 if length(MRCA) > 1
#                     MRCA = intersect(MRCA...)
#                     MRCA = MRCA[MRCA .!= [0]]
#                     MRCA = maximum(MRCA)
#                 else
#                     MRCA = lineage
#                 end
#                 (t_creation,) = execute(dbTemp,"SELECT t_creation FROM $(s)strains
#                                         WHERE $(s)strain_id = $(MRCA)")
#                 execute(dbOutput, "BEGIN TRANSACTION")
#                 for daughter in daughters
#                     execute(dbOutput, "INSERT INTO MRCA_$(s)strains VALUES (?,?,?,?,?)",
#                                             (time,daughter,lineage,MRCA,t_creation.t_creation))
#                 end
#                 strainpairs = collect(combinations(strain_idx,2))
#                 for pair in strainpairs
#                     strainID1 = pair[1]
#                     strainID2 = pair[2]
#                     daughters = strains[pair]
#                     MRCA = ancestors[[strainID1,strainID2]]
#                     if length(MRCA) > 1
#                         MRCA = intersect(MRCA...)
#                         MRCA = MRCA[MRCA .!= [0]]
#                         MRCA = maximum(MRCA)
#                     else
#                         MRCA = lineage
#                     end
#                     (t_creation,) = execute(dbTemp,"SELECT t_creation FROM $(s)strains
#                                             WHERE $(s)strain_id = $(MRCA)")
#                     execute(dbOutput, "INSERT INTO pairwise_MRCA_$(s)strains VALUES (?,?,?,?,?,?)",
#                             (time,daughters[1],daughters[2],
#                             lineage,MRCA,t_creation.t_creation))
#                 end
#                 execute(dbOutput, "COMMIT")
#             end
#         end
#     end
# end


keepStrains, keepVstrains = keepWhichStrains(bthreshold, vthreshold)
orderTreeStrains("microbe", keepStrains, keepVstrains)
orderTreeStrains("virus", keepStrains, keepVstrains)
findExtinctionTimes("virus", keepVstrains)
findExtinctionTimes("microbe", keepStrains)



# function createindices()
#     println("(Creating run_id indices...)")
#     db = SQLite.DB(dbOutputPath)
#     table_names = [table_name for (table_name,) in execute(
#         db, "SELECT name FROM sqlite_schema
#         WHERE type='table' ORDER BY name;")]
#     # execute(db, "BEGIN TRANSACTION")
#     dbOutput = SQLite.DB(dbOutputPath)
#     for table_name in table_names
#         # cols = [info.name for info in execute(db,"PRAGMA table_info($(table_name))")]
#         if in(table_name,["tree_bstrain_order"])
#             execute(dbOutput, "CREATE INDEX $(table_name)_index ON $(table_name) (bstrain_id)")
#         end
#         if in(table_name,["tree_vstrain_order"])
#             execute(dbOutput, "CREATE INDEX $(table_name)_index ON $(table_name) (vstrain_id)")
#         end
#         if in(table_name,["tree_babundance","tree_vabundance"])
#             execute(dbOutput, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
#         end
#     end
#     # execute(dbSim, "COMMIT")
# end
# createindices()
if isfile(joinpath(SCRIPT_PATH, "$(dbSimDir)", "runID-cID$(combo_id).txt"))
    rm(joinpath(SCRIPT_PATH, "$(dbSimDir)", "runID-cID$(combo_id).txt"))
end
open(joinpath(SCRIPT_PATH, "$(dbSimDir)", "runID-cID$(combo_id).txt"), "w") do f
    println(f, "$(run_id)")
end
println("Complete!")
