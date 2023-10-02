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
dbOutputPath = joinpath("mrca_output.sqlite") # cluster

dbOutput = SQLite.DB(dbOutputPath)

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
execute(dbOutput, "CREATE TABLE expected_MRCA_bstrains 
                    (t REAL , lineage INTEGER, exp_t_creation REAL)")
execute(dbOutput, "CREATE TABLE expected_MRCA_vstrains 
                    (t REAL , lineage INTEGER, exp_t_creation REAL)")                   
execute(dbOutput, "CREATE TABLE expected_pairwise_MRCA_bstrains 
                    (t REAL, lineage INTEGER, exp_t_creation REAL)")
execute(dbOutput, "CREATE TABLE expected_pairwise_MRCA_vstrains 
                    (t REAL, lineage INTEGER, exp_t_creation REAL)")



# Create temporary database that is a copy of the main database at the run_id value of the script's argument
dbTemp = SQLite.DB()
#dbTemp = SQLite.DB("/Volumes/Yadgah/test.sqlite") # local
# dbTemp = SQLite.DB(dbSimPath) # local
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE bstrains (t_creation REAL, bstrain_id INTEGER, parent_bstrain_id INTEGER, infecting_vstrain_id INTEGER)")
execute(dbTemp, "CREATE TABLE vstrains (t_creation REAL, vstrain_id INTEGER, parent_vstrain_id INTEGER, infected_bstrain_id INTEGER)")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp,"ATTACH DATABASE '$(dbSimPath)' as dbSim")

execute(dbTemp,"INSERT INTO bstrains (t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id)
SELECT t_creation, bstrain_id,parent_bstrain_id,infecting_vstrain_id
FROM dbSim.bstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vstrains (t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id)
SELECT t_creation, vstrain_id,parent_vstrain_id,infected_bstrain_id
FROM dbSim.vstrains WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO babundance (t, bstrain_id, abundance)
SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp,"INSERT INTO vabundance (t, vstrain_id, abundance)
SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")


execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX babundance_index ON babundance (t,bstrain_id)")
execute(dbTemp, "CREATE INDEX vabundance_index ON vabundance (t,vstrain_id)")
execute(dbTemp, "CREATE INDEX bstrains_index ON bstrains (bstrain_id,parent_bstrain_id)")
execute(dbTemp, "CREATE INDEX vstrains_index ON vstrains (vstrain_id,parent_vstrain_id)")
execute(dbTemp, "COMMIT")


function orderTreeStrains(strainType::String)
    if strainType == "virus"
        s = "v"
        initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)strains WHERE parent_$(s)strain_id = 0 ORDER BY $(s)strain_id")]
        # strainTreeOrder = DataFrame(bstrain_id = Int64[0,1], tree_bstrain_id = Int64[0,1])
        lineages = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
            ORDER BY t_creation")
            lineages[strain_id] = searchLineage(strain_id,strainType)
        end
        MRCAtimeSeries(lineages,strainType)
    end
    if strainType == "microbe"
        s = "b"
        initialStrains =[strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
        FROM $(s)strains WHERE parent_$(s)strain_id = 0 ORDER BY $(s)strain_id")]
        # strainTreeOrder = DataFrame(bstrain_id = Int64[0,1], tree_bstrain_id = Int64[0,1])
        lineages = Dict{Int64,Array{Int64,1}}()
        for (strain_id,) in execute(dbTemp,"SELECT DISTINCT $(s)strain_id FROM $(s)strains
            ORDER BY t_creation")
            lineages[strain_id] = searchLineage(strain_id,strainType)
        end
        MRCAtimeSeries(lineages,strainType)
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
        pairTMRCA = 0
        TMRCA = 0
        nTMRCA = 0
        nPairTMRCA = 0
        for (time,) in execute(dbTemp,"SELECT DISTINCT t FROM $(s)abundance ORDER BY t")
            strains = [strain_id for (strain_id,) in execute(dbTemp,"SELECT $(s)strain_id
            FROM $(s)abundance WHERE t = $(time) ORDER BY $(s)strain_id")]
            println("$(strainType) MRCA at time $(time)")
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
                    TMRCA += t_creation.t_creation
                    nTMRCA += 1
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
                    pairTMRCA += t_creation.t_creation
                    nPairTMRCA += 1
                end
                execute(dbOutput, "INSERT INTO expected_MRCA_$(s)strains VALUES (?,?,?)",
                    (time, lineage, TMRCA / nTMRCA))
                TMRCA = 0
                nTMRCA = 0
                execute(dbOutput, "INSERT INTO expected_pairwise_MRCA_$(s)strains VALUES (?,?,?)",
                    (time, lineage, pairTMRCA / nPairTMRCA))
                pairTMRCA = 0 
                nPairTMRCA = 0
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
# newAbundanceTimeSeries("microbe")
orderTreeStrains("virus")
# newAbundanceTimeSeries("virus")
# findExtinctionTimes("microbe")
# findExtinctionTimes("virus")


# function createindices()
#     println("(Creating run_id indices...)")
#     db = SQLite.DB(dbOutputPath)
#     execute(db, "BEGIN TRANSACTION")
#     for (table_name,) in execute(
#         db, "SELECT name FROM sqlite_schema
#         WHERE type='table' ORDER BY name;")
#         # cols = [info.name for info in execute(db,"PRAGMA table_info($(table_name))")]
#         if in(table_name,["tree_bstrain_order"])
#             execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (bstrain_id)")
#         end
#         if in(table_name,["tree_vstrain_order"])
#             execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (vstrain_id)")
#         end
#         if in(table_name,["tree_babundance","tree_vabundance"])
#             execute(db, "CREATE INDEX $(table_name)_index ON $(table_name) (t)")
#         end
#     end
#     execute(db, "COMMIT")
# end
# createindices()


println("Complete!")
