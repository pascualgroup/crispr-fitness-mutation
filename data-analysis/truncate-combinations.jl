#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH, "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("db-truncate-combinations.sqlite") # cluster
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

# dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/21_MOI3", "sweep_db_gathered.sqlite")
# dbOutputPath = joinpath("/Volumes/Yadgah", "truncate-test.sqlite")

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE summary (t REAL, viral_abundance INTEGER,  microbial_abundance INTEGER, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE bgrowthrates (bstrain_id INTEGER,  growth_rate REAL, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER,  abundance INTEGER, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE bextinctions (t_extinction REAL, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE bspacers (bstrain_id INTEGER,  spacer_id INTEGER, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE bstrains (t_creation REAL, bstrain_id INTEGER,  parent_bstrain_id INTEGER, infecting_vstrain_id INTEGER, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER,  abundance REAL, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE vextinctions (t_extinction REAL, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE vpspacers (vstrain_id INTEGER,  spacer_id INTEGER, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE vstrains (t_creation REAL, vstrain_id INTEGER,  parent_vstrain_id INTEGER, infected_bstrain_id INTEGER, run_id INTEGER)")

cIDs = [cID for (cID,) 
            in execute(dbSim,"SELECT combo_id FROM param_combos 
                    WHERE microbe_carrying_capacity = 400000
                    AND init_bcomm_function in (1,2)
                    AND microbe_mutation_prob in (0) 
                    AND evofunctionScale in (0,3)")]
append!(cIDs,[cID for (cID,) 
            in execute(dbSim,"SELECT combo_id FROM param_combos 
                    WHERE microbe_carrying_capacity = 400000
                    AND init_bcomm_function in (1,2) 
                    AND microbe_mutation_prob in (1) 
                    AND evofunctionScale in (3)")])
runIDs = [runID for (runID,) 
            in execute(dbSim, "SELECT run_id FROM runs 
                        WHERE combo_id in ($(join(cIDs, ", ")))")]

execute(dbOutput, "ATTACH DATABASE '$(dbSimPath)' as dbSim")

execute(dbOutput, "BEGIN TRANSACTION")

execute(dbOutput, "INSERT INTO summary (t, viral_abundance, microbial_abundance, run_id) 
                    SELECT t, viral_abundance, microbial_abundance, run_id 
                    FROM dbSim.summary WHERE run_id in ($(join(runIDs, ", ")));")
println("...summary complete...")

execute(dbOutput, "INSERT INTO bgrowthrates (bstrain_id, growth_rate, run_id) 
                    SELECT bstrain_id, growth_rate, run_id 
                    FROM dbSim.bgrowthrates WHERE run_id in ($(join(runIDs, ", ")));")
println("...bgrowthrates complete...")

execute(dbOutput, "INSERT INTO babundance (t, bstrain_id, abundance, run_id) 
                    SELECT t, bstrain_id, abundance, run_id 
                    FROM dbSim.babundance WHERE run_id in ($(join(runIDs, ", ")));")
println("...babundance complete...")

execute(dbOutput, "INSERT INTO bspacers (bstrain_id, spacer_id, run_id) 
                    SELECT bstrain_id, spacer_id, run_id 
                    FROM dbSim.bspacers WHERE run_id in ($(join(runIDs, ", ")));")
println("...bspacers complete...")

execute(dbOutput, "INSERT INTO bstrains (t_creation, bstrain_id, parent_bstrain_id, infecting_vstrain_id, run_id) 
                    SELECT t_creation, bstrain_id, parent_bstrain_id, infecting_vstrain_id, run_id
                    FROM dbSim.bstrains WHERE run_id in ($(join(runIDs, ", ")));")
println("...bstrains complete...")

execute(dbOutput, "INSERT INTO vabundance (t, vstrain_id, abundance, run_id) 
                    SELECT t, vstrain_id, abundance, run_id 
                    FROM dbSim.vabundance WHERE run_id in ($(join(runIDs, ", ")));")
println("...vabundance complete...")

execute(dbOutput, "INSERT INTO vpspacers (vstrain_id, spacer_id, run_id) 
                    SELECT vstrain_id, spacer_id, run_id 
                    FROM dbSim.vpspacers WHERE run_id in ($(join(runIDs, ", ")));")
println("...vpspacers complete...")

execute(dbOutput, "INSERT INTO vstrains (t_creation, vstrain_id, parent_vstrain_id, infected_bstrain_id, run_id) 
                    SELECT t_creation, vstrain_id, parent_vstrain_id, infected_bstrain_id, run_id
                    FROM dbSim.vstrains WHERE run_id in ($(join(runIDs, ", ")));")
println("...vstrains complete...")

execute(dbOutput, "COMMIT")

execute(dbOutput, "DETACH DATABASE dbSim")

DataFrame(execute(dbSim, "SELECT * FROM param_combos")) |> SQLite.load!(dbOutput, "param_combos", ifnotexists=true)
DataFrame(execute(dbSim, "SELECT * FROM runs")) |> SQLite.load!(dbOutput, "runs", ifnotexists=true)


execute(dbOutput, "BEGIN TRANSACTION")

for (run_id,) in execute(dbSim, "SELECT run_id FROM runs WHERE run_id in ($(join(runIDs, ", "))) ORDER BY run_id")
    println(run_id)
    DF = DataFrame(execute(
        dbSim,
        "SELECT t, viral_abundance FROM summary WHERE run_id = $(run_id) AND viral_abundance = 0 ORDER by t DESC"
    ))
    if size(DF)[1] != 0
        (tExt,) = execute(
            dbSim,
            "SELECT t_extinction FROM vextinctions 
            WHERE run_id = $(run_id) ORDER BY t_extinction DESC LIMIT 1")
        execute(dbOutput, "INSERT INTO vextinctions VALUES ($(tExt[1]),$(run_id))")
    end
    DF = DataFrame(execute(
        dbSim,
        "SELECT t, microbial_abundance FROM summary WHERE run_id = $(run_id) AND microbial_abundance = 0 ORDER by t DESC"
    ))
    if size(DF)[1] != 0
        (tExt,) = execute(
            dbSim,
            "SELECT t_extinction FROM bextinctions 
            WHERE run_id = $(run_id) ORDER BY t_extinction DESC LIMIT 1")
        execute(dbOutput, "INSERT INTO bextinctions VALUES ($(tExt[1]),$(run_id))")
    end
end
execute(dbOutput, "COMMIT")

println("Complete!")

