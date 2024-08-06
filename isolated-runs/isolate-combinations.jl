#!/usr/bin/env julia
println("(Julia compilation delay...)")

using SQLite
import SQLite.DBInterface.execute

## Define Paths ##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
dbSimPath = joinpath(SCRIPT_PATH,"..","simulation","sweep_db_gathered.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah/crispr-sweep-13-1-2021/simulation/sweep_db_gathered.sqlite") # local
dbSim = SQLite.DB(dbSimPath)
##

function main()
    for (cID,) in execute(dbSim, "SELECT combo_id FROM param_combos ORDER BY combo_id")
        println("Isolating combination $(cID)...")
        run_ids = [runID for (runID,) 
                    in execute(dbSim, 
                    "SELECT run_id FROM runs 
                    WHERE combo_id = $(cID)
                    ORDER BY combo_id")]
            if !ispath(joinpath(SCRIPT_PATH,"combo-isolates"))
                mkpath(joinpath(SCRIPT_PATH,"combo-isolates"))
            end
        dbOutput = SQLite.DB(joinpath(SCRIPT_PATH, "combo-isolates", "comboID-$(cID).sqlite"))
        println("...attaching database...")
        execute(dbOutput, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
        println("Database Attached")

        for table_name in ["babundance","bspacers","bgrowthrates","bstrains","summary","vabundance","vpspacers","vstrains","bextinctions","vextinctions"]
            if length([(name,) for name in execute(dbSim, "SELECT name FROM sqlite_master WHERE type='table' AND name='$(table_name)';")]) > 0
                println("Table Name: $(table_name)")
                tableCols = ["$(table_info.name)" for table_info in execute(dbSim,"PRAGMA table_info($(table_name))")]
                tableColsType = ["$(table_info.name) $(table_info.type)" for table_info in execute(dbSim,"PRAGMA table_info($(table_name))")]
                numCols = length(tableCols)
                colStmt = join(tableCols,", ")
                colTypeStmt = join(tableColsType,", ")
                execute(dbOutput, "BEGIN TRANSACTION")
                println("Creating Table")
                execute(dbOutput, "CREATE TABLE $(table_name) ($(colTypeStmt...))")
                println("Loading $(table_name) data")
                execute(dbOutput, "INSERT INTO $(table_name)($(colStmt)) SELECT $(colStmt)
                FROM dbSim.$(table_name) WHERE run_id in ($(join(run_ids,", ")));")
                println("Table Created")
                execute(dbOutput, "COMMIT")
            end
        end

        println("Table Name: param_combos & runs")
        tableNamesTypes = ["$(table_info.name) $(table_info.type)" for table_info in execute(dbSim,"PRAGMA table_info(param_combos)")]
        tableNamesTypes = join(tableNamesTypes,", ")
        println("...creating tables...")
        execute(dbOutput, "CREATE TABLE runs (run_id INTEGER, combo_id INTEGER, replicate INTEGER, rng_seed INTEGER)")
        execute(dbOutput, "CREATE TABLE param_combos ($(tableNamesTypes...))")
        println("Tables Created")
        tableNames = ["$(table_info.name)" for table_info in execute(dbSim,"PRAGMA table_info(param_combos)")]
        tableNames = join(tableNames,", ")
        execute(dbOutput, "BEGIN TRANSACTION")
        println("...loading param_combos & runs data...")
        execute(dbOutput,"INSERT INTO param_combos($(tableNames)) SELECT * FROM dbSim.param_combos")
        execute(dbOutput,"INSERT INTO runs (run_id, combo_id, replicate) SELECT run_id, combo_id, replicate FROM dbSim.runs")
        execute(dbOutput, "COMMIT")
    end
end

main()
println("Combinations successfully isolated!")
