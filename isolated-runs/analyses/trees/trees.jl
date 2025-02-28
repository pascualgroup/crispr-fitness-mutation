#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute
using Combinatorics


run_id = ARGS[1]

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH,"..","..","..","simulation","sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("trees_output.sqlite") # cluster

# dbSimPath = joinpath("/Volumes/Yadgah/sweep_db_gathered.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah","run_id1455_combo73_replicate15.sqlite") # local
# dbSimPath = joinpath("/Users/armun/Dropbox/Current/Projects/microbe-virus-crispr/stochastic-crispr",
#                             "runID146-c5-r6.sqlite") # local
# dbSimPath = joinpath("/Volumes/Yadgah/crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite")
# dbSimPath = joinpath("/Volumes/Yadgah", "sylvain-martin-collab/12_MOI3/isolates/runID4209-c15-r9/runID4209-c15-r9.sqlite") # local 
# dbOutputPath = joinpath("/Volumes/Yadgah/trees_output.sqlite") # local
# dbOutputPath = joinpath("/Users/armun/Dropbox/Current/Projects/microbe-virus-crispr/stochastic-crispr/trees_output.sqlite") # local
##

dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE tree_babundance (t REAL, tree_bstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE tree_vabundance (t REAL, tree_vstrain_id INTEGER, abundance INTEGER)")
execute(dbOutput, "CREATE TABLE bstrain_creation_extinction (bstrain_id INTEGER, t_creation REAL, t_extinction REAL)")
execute(dbOutput, "CREATE TABLE tree_bstrain_creation_extinction
(tree_bstrain_id INTEGER, t_creation REAL, t_extinction REAL,
tree_parent_bstrain_id INTEGER, tree_infecting_vstrain_id)")
execute(dbOutput, "CREATE TABLE vstrain_creation_extinction (vstrain_id INTEGER, t_creation REAL, t_extinction REAL)")
execute(dbOutput, "CREATE TABLE tree_vstrain_creation_extinction
(tree_vstrain_id INTEGER, t_creation REAL, t_extinction REAL,
tree_parent_vstrain_id INTEGER, tree_infected_bstrain_id INTEGER)")
execute(dbOutput, "CREATE TABLE MRCA_bstrains
                    (t REAL, bstrain_id INTEGER, lineage INTEGER, MRCA_bstrain_id INTEGER,
                        t_creation REAL)")
execute(dbOutput, "CREATE TABLE MRCA_vstrains
                    (t REAL, vstrain_id INTEGER, lineage INTEGER, MRCA_vstrain_id INTEGER,
                        t_creation REAL)")
execute(dbOutput, "CREATE TABLE pairwise_MRCA_bstrains
                    (t REAL, bstrain_id_sample_1 INTEGER, bstrain_id_sample_2 INTEGER,
                    lineage INTEGER, MRCA_bstrain_id INTEGER, t_creation REAL)")
execute(dbOutput, "CREATE TABLE pairwise_MRCA_vstrains
                    (t REAL, vstrain_id_sample_1 INTEGER, vstrain_id_sample_2 INTEGER,
                    lineage INTEGER, MRCA_vstrain_id INTEGER, t_creation REAL)")



# Create temporary database that is a copy of the main database at the run_id value of the script's argument
dbTemp = SQLite.DB()
#dbTemp = SQLite.DB("/Volumes/Yadgah/test.sqlite") # local
# dbTemp = SQLite.DB(dbSimPath) # local
# execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
# execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE bstrains (t_creation REAL, bstrain_id INTEGER, parent_bstrain_id INTEGER, infecting_vstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE vstrains (t_creation REAL, vstrain_id INTEGER, parent_vstrain_id INTEGER, infected_bstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE bgrowthrates (bstrain_id INTEGER, growth_rate REAL)")
execute(dbTemp, "CREATE TABLE bextinctions (bstrain_id INTEGER, t_extinction REAL)")
execute(dbTemp, "CREATE TABLE vextinctions (vstrain_id INTEGER, t_extinction REAL)")
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")
# execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
# execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bstrains (t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id)
SELECT t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id
FROM dbSim.bstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vstrains (t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id)
SELECT t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id
FROM dbSim.vstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bgrowthrates (bstrain_id, growth_rate)
SELECT bstrain_id, growth_rate
FROM dbSim.bgrowthrates WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO bextinctions (bstrain_id, t_extinction)
SELECT bstrain_id, t_extinction
FROM dbSim.bextinctions WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vextinctions (vstrain_id, t_extinction)
SELECT vstrain_id, t_extinction
FROM dbSim.vextinctions WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance)
SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance)
SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX babundance_index ON babundance (t,bstrain_id,abundance)")
execute(dbTemp, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id, abundance)")
execute(dbTemp, "CREATE INDEX bstrains_index ON bstrains (t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id)")
execute(dbTemp, "CREATE INDEX vstrains_index ON vstrains (t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id)")
execute(dbTemp, "CREATE INDEX bgrowth_index ON bgrowthrates (bstrain_id)")
execute(dbTemp, "COMMIT")


function orderTreeStrains(strainType::String)
    if strainType == "virus"
        s = "v"
        initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)strains WHERE parent_$(s)strain_id = 0 ORDER BY $(s)strain_id")]
        strainTreeOrder = DataFrame(bstrain_id = Int64[0,1], tree_bstrain_id = Int64[0,1])
        lineages = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
            ORDER BY t_creation")
            lineages[strain_id] = searchLineage(strain_id,strainType)
        end
        initialtreestrainIDs = Int64[0, 1]
        for initialStrain in initialStrains
            if initialStrain == 1
                continue
            end
            newTreeID = initialtreestrainIDs[initialStrain] + length(lineages[initialStrain-1])
            push!(initialtreestrainIDs,newTreeID)
        end
        strainTreeOrder = DataFrame(vstrain_id = Int64[0,initialStrains...],
                                        tree_vstrain_id = initialtreestrainIDs)
        for initialStrain in initialStrains
            MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
            FROM $(s)strains WHERE parent_$(s)strain_id = 1 ORDER BY $(s)strain_id")]
            while length(MRCA_ids) > 0
                println("Finding positions for MRCA: $(MRCA_ids)")
                nextGeneration!(strainTreeOrder,MRCA_ids,lineages,strainType)
                MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
                FROM $(s)strains WHERE parent_$(s)strain_id in ($(join(MRCA_ids,", ")))
                ORDER BY $(s)strain_id")]
            end
        end
        strainTreeOrder |> SQLite.load!(dbOutput,"tree_vstrain_order",ifnotexists=true)
        MRCAtimeSeries(lineages,strainType)
    end
    if strainType == "microbe"
        s = "b"
        initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)strains WHERE parent_$(s)strain_id = 0 ORDER BY $(s)strain_id")]
        initialStrainsGrowthOrdered =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)growthrates WHERE $(s)strain_id in ($(join(initialStrains,", "))) ORDER BY growth_rate DESC")]
        strainTreeOrder = DataFrame(bstrain_id = Int64[0,1], tree_bstrain_id = Int64[0,1])
        lineages = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
            ORDER BY t_creation")
            lineages[strain_id] = searchLineage(strain_id,strainType)
        end
        initialtreestrainIDs = Int64[0, 1]
        # for initialStrain in initialStrains
        for strainIDX in collect(2:1:length(initialStrainsGrowthOrdered))
            # if initialStrain == 1
            #     continue
            # end
            strainBefore = initialStrainsGrowthOrdered[strainIDX-1]
            newTreeID = initialtreestrainIDs[strainIDX] + length(lineages[strainBefore])
            push!(initialtreestrainIDs,newTreeID)
        end
        strainTreeOrder = DataFrame(bstrain_id = Int64[0,initialStrainsGrowthOrdered...],
                                        tree_bstrain_id = initialtreestrainIDs)
        for initialStrain in initialStrainsGrowthOrdered
            MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
            FROM $(s)strains WHERE parent_$(s)strain_id = $(initialStrain) ORDER BY $(s)strain_id")]
            while length(MRCA_ids) > 0
                println("Finding positions for MRCA: $(MRCA_ids)")
                nextGeneration!(strainTreeOrder,MRCA_ids,lineages,strainType)
                MRCA_ids = [MRCA_id for (MRCA_id,) in execute(dbTemp,"SELECT $(s)strain_id
                FROM $(s)strains WHERE parent_$(s)strain_id in ($(join(MRCA_ids,", ")))
                ORDER BY $(s)strain_id")]
            end
        end
        strainTreeOrder |> SQLite.load!(dbOutput,"tree_bstrain_order",ifnotexists=true)
        MRCAtimeSeries(lineages,strainType)
    end
end


function nextGeneration!(strainTreeOrder::DataFrame,
    MRCA_ids::Array{Int64,1},lineages::Dict,strainType::String)
    # MRCA_ids = sort([Int64(key) for key in keys(lineages)])
    if strainType == "virus"
        s = "v"
        parentStrainIDs = Int64[]
        for MRCA_id in MRCA_ids
            (strain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(MRCA_id) ORDER BY parent_$(s)strain_id")
            union!(parentStrainIDs,strain.parent_vstrain_id)
        end
        sort!(parentStrainIDs)
        # println("Parents of MRCA are: $(parentStrainIDs)")
        treePositions!(strainTreeOrder,parentStrainIDs,lineages,strainType)
    end
    if strainType == "microbe"
        s = "b"
        parentStrainIDs = Int64[]
        for MRCA_id in MRCA_ids
            (strain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(MRCA_id) ORDER BY parent_$(s)strain_id")
            union!(parentStrainIDs,strain.parent_bstrain_id)
        end
        sort!(parentStrainIDs)
        # println("Parents of MRCA are: $(parentStrainIDs)")
        treePositions!(strainTreeOrder,parentStrainIDs,lineages,strainType)
    end
end

function treePositions!(strainTreeOrder::DataFrame,
    parentStrainIDs::Array{Int64,1},lineages::Dict,strainType::String)
    if strainType == "virus"
        s = "v"
        for parentID in parentStrainIDs
            # println("Locating descendants of strain $(parentID)")
            subMRCA_ids = [strain for (strain,) in execute(dbTemp, "SELECT $(s)strain_id FROM $(s)strains
            WHERE parent_$(s)strain_id = $(parentID) ORDER BY $(s)strain_id")]
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
            WHERE parent_$(s)strain_id = $(parentID) ORDER BY $(s)strain_id")]
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

function newAbundanceTimeSeries(strainType::String)
    if strainType == "virus"
        s = "v"
        for (time,) in execute(dbTemp,"SELECT DISTINCT t
            FROM $(s)abundance ORDER BY t")
            println("Changing strain IDs at time = $(time)")
            for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)abundance
                WHERE t = $(time) ORDER BY $(s)strain_id")
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
            println("Changing strain IDs at time = $(time)")
            for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)abundance
                WHERE t = $(time) ORDER BY $(s)strain_id")
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

function findExtinctionTimes(strainType::String)
    if strainType == "virus"
        s = "v"
        for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)extinctions
                                    ORDER BY $(s)strain_id")
            (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            (tExt,) = execute(dbTemp, "SELECT t_extinction FROM $(s)extinctions
            WHERE $(s)strain_id = $(strain)")
            execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
            VALUES (?,?,?)",(strain,t.t_creation,tExt.t_extinction))
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
            VALUES (?,?,?,?,?)",(tree.tree_vstrain_id,t.t_creation,tExt.t_extinction,
            pstrain.tree_vstrain_id,ibstrain.tree_bstrain_id))
        end
    end
    if strainType == "microbe"
        s = "b"
        println("Searching for strain extinctions")
        for (strain,) in execute(dbTemp,"SELECT $(s)strain_id FROM $(s)extinctions 
                                    ORDER BY $(s)strain_id")
            (t,) = execute(dbTemp, "SELECT t_creation FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            (tExt,) = execute(dbTemp, "SELECT t_extinction FROM $(s)extinctions
            WHERE $(s)strain_id = $(strain)")
            execute(dbOutput,"INSERT INTO $(s)strain_creation_extinction
            VALUES (?,?,?)",(strain,t.t_creation,tExt.t_extinction))
            (pstrain,) = execute(dbTemp, "SELECT parent_$(s)strain_id FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            (ivstrain,) = execute(dbTemp, "SELECT infecting_vstrain_id FROM $(s)strains
            WHERE $(s)strain_id = $(strain)")
            (pstrain,) = execute(dbOutput, "SELECT tree_$(s)strain_id
            FROM tree_$(s)strain_order WHERE $(s)strain_id = $(pstrain.parent_bstrain_id)")
            (ivstrain,) = execute(dbOutput, "SELECT tree_vstrain_id
            FROM tree_vstrain_order WHERE vstrain_id = $(ivstrain.infecting_vstrain_id)")
            (tree,) = execute(dbOutput, "SELECT tree_$(s)strain_id
            FROM tree_$(s)strain_order WHERE $(s)strain_id = $(strain)")
            execute(dbOutput,"INSERT INTO tree_$(s)strain_creation_extinction
            VALUES (?,?,?,?,?)",(tree.tree_bstrain_id,t.t_creation,tExt.t_extinction,
            pstrain.tree_bstrain_id,ivstrain.tree_vstrain_id))
        end
    end
end

function searchLineage(MRCA_id::Int64,strainType::String)
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
            WHERE parent_$(type)strain_id in ($(id)) ORDER BY $(type)strain_id")]
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

function searchAncestry(daughter_id::Int64,strainType::String)
    if strainType == "virus"
        type = "v"
    end
    if strainType == "microbe"
        type = "b"
    end
    ancestors = Int64[daughter_id]
    ancestor_id = daughter_id
    while ancestor_id > 0
        ancestor_id = [id for (id,) in
        execute(dbTemp,"SELECT parent_$(type)strain_id FROM $(type)strains
        WHERE $(type)strain_id in ($(ancestor_id)) ORDER BY $(type)strain_id")][1]
        push!(ancestors,ancestor_id)
    end
    return sort(ancestors)
end

function MRCAtimeSeries(lineages::Dict{Int64,Array{Int64,1}},strainType::String)
    if strainType == "virus"
        s = "v"
        ancestry = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
            ORDER BY t_creation")
            anc = searchAncestry(strain_id,strainType)
            ancestry[strain_id] = anc[anc.!=strain_id]
        end
        initialStrains = [strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)strains WHERE parent_$(s)strain_id = 0 ORDER BY $(s)strain_id")]
        for (time,) in execute(dbTemp,"SELECT DISTINCT t FROM $(s)abundance ORDER BY t")
            strains = [strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
            FROM $(s)abundance WHERE t = $(time) ORDER BY $(s)strain_id")]
            if length(strains) == 0
                continue
            end
            ancestors = map(x->ancestry[x],strains)
            for lineage in initialStrains
                strain_idx = findall(x-> lineage in x, ancestors)
                daughters = strains[strain_idx]
                MRCA = ancestors[strain_idx]
                if length(MRCA) > 1
                    MRCA = intersect(MRCA...)
                    MRCA = MRCA[MRCA .!= [0]]
                    MRCA = maximum(MRCA)
                else
                    MRCA = lineage
                end
                (t_creation,) = execute(dbTemp,"SELECT t_creation FROM $(s)strains
                                        WHERE $(s)strain_id = $(MRCA)")
                execute(dbOutput, "BEGIN TRANSACTION")
                for daughter in daughters
                    execute(dbOutput, "INSERT INTO MRCA_$(s)strains VALUES (?,?,?,?,?)",
                                            (time,daughter,lineage,MRCA,t_creation.t_creation))
                end
                strainpairs = collect(combinations(strain_idx,2))
                for pair in strainpairs
                    strainID1 = pair[1]
                    strainID2 = pair[2]
                    daughters = strains[pair]
                    MRCA = ancestors[[strainID1,strainID2]]
                    if length(MRCA) > 1
                        MRCA = intersect(MRCA...)
                        MRCA = MRCA[MRCA .!= [0]]
                        MRCA = maximum(MRCA)
                    else
                        MRCA = lineage
                    end
                    (t_creation,) = execute(dbTemp,"SELECT t_creation FROM $(s)strains
                                            WHERE $(s)strain_id = $(MRCA)")
                    execute(dbOutput, "INSERT INTO pairwise_MRCA_$(s)strains VALUES (?,?,?,?,?,?)",
                            (time,daughters[1],daughters[2],
                            lineage,MRCA,t_creation.t_creation))
                end
                execute(dbOutput, "COMMIT")
            end
        end
    end
    if strainType == "microbe"
        s = "b"
        ancestry = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
            ORDER BY t_creation")
            anc = searchAncestry(strain_id,strainType)
            ancestry[strain_id] = anc[anc.!=strain_id]
        end
        initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)strains WHERE parent_$(s)strain_id = 0 ORDER BY $(s)strain_id")]
        for (time,) in execute(dbTemp,"SELECT DISTINCT t FROM $(s)abundance ORDER BY t")
            strains = [strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
            FROM $(s)abundance WHERE t = $(time) ORDER BY $(s)strain_id")]
            if length(strains) == 0
                continue
            end
            ancestors = map(x->ancestry[x],strains)
            for lineage in initialStrains
                strain_idx = findall(x-> lineage in x, ancestors)
                daughters = strains[strain_idx]
                MRCA = ancestors[strain_idx]
                if length(MRCA) > 1
                    MRCA = intersect(MRCA...)
                    MRCA = MRCA[MRCA .!= [0]]
                    MRCA = maximum(MRCA)
                else
                    MRCA = lineage
                end
                (t_creation,) = execute(dbTemp,"SELECT t_creation FROM $(s)strains
                                        WHERE $(s)strain_id = $(MRCA)")
                execute(dbOutput, "BEGIN TRANSACTION")
                for daughter in daughters
                    execute(dbOutput, "INSERT INTO MRCA_$(s)strains VALUES (?,?,?,?,?)",
                                            (time,daughter,lineage,MRCA,t_creation.t_creation))
                end
                strainpairs = collect(combinations(strain_idx,2))
                for pair in strainpairs
                    strainID1 = pair[1]
                    strainID2 = pair[2]
                    daughters = strains[pair]
                    MRCA = ancestors[[strainID1,strainID2]]
                    if length(MRCA) > 1
                        MRCA = intersect(MRCA...)
                        MRCA = MRCA[MRCA .!= [0]]
                        MRCA = maximum(MRCA)
                    else
                        MRCA = lineage
                    end
                    (t_creation,) = execute(dbTemp,"SELECT t_creation FROM $(s)strains
                                            WHERE $(s)strain_id = $(MRCA)")
                    execute(dbOutput, "INSERT INTO pairwise_MRCA_$(s)strains VALUES (?,?,?,?,?,?)",
                            (time,daughters[1],daughters[2],
                            lineage,MRCA,t_creation.t_creation))
                end
                execute(dbOutput, "COMMIT")
            end
        end
    end
end

orderTreeStrains("microbe")
newAbundanceTimeSeries("microbe")
orderTreeStrains("virus")
newAbundanceTimeSeries("virus")
findExtinctionTimes("microbe")
findExtinctionTimes("virus")


function createindices()
    println("(Creating run_id indices...)")
    db = SQLite.DB(dbOutputPath)
    table_names = [table_name for (table_name,) in execute(
        db, "SELECT name FROM sqlite_schema
        WHERE type='table' ORDER BY name;")]
    # execute(db, "BEGIN TRANSACTION")
    dbOutput = SQLite.DB(dbOutputPath)
    for table_name in table_names
        # cols = [info.name for info in execute(db,"PRAGMA table_info($(table_name))")]
        if in(table_name,["tree_bstrain_order"])
            execute(dbOutput, "CREATE INDEX $(table_name)_index ON $(table_name) (bstrain_id)")
        end
        if in(table_name,["tree_vstrain_order"])
            execute(dbOutput, "CREATE INDEX $(table_name)_index ON $(table_name) (vstrain_id)")
        end
        if in(table_name,["tree_babundance","tree_vabundance"])
            execute(dbOutput, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
        end
    end
    # execute(dbSim, "COMMIT")
end
createindices()


println("Complete!")
