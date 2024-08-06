#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH, "..", "simulation", "sweep_db_gathered.sqlite") # cluster
dbOutputPath = joinpath("db-truncate.sqlite") # cluster
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE bgrowthrates (bstrain_id INTEGER,  growth_rate REAL, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE summary (t INTEGER, viral_abundance INTEGER,  microbial_abundance INTEGER, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE babundance (t INTEGER, bstrain_id INTEGER,  abundance INTEGER, run_id INTEGER)")
execute(dbOutput, "CREATE TABLE vextinctions (t_extinction REAL, run_id INTEGER)")


execute(dbOutput, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbOutput, "BEGIN TRANSACTION")
execute(dbOutput, "INSERT INTO bgrowthrates (bstrain_id, growth_rate, run_id) SELECT bstrain_id, growth_rate, run_id FROM dbSim.bgrowthrates;")
println("...bgrowthrates complete...")
execute(dbOutput, "INSERT INTO summary (t, viral_abundance, microbial_abundance, run_id) SELECT t, viral_abundance, microbial_abundance, run_id 
                    FROM dbSim.summary;")
println("...summary complete...")
execute(dbOutput, "INSERT INTO babundance (t, bstrain_id, abundance, run_id) SELECT t, bstrain_id, abundance, run_id FROM dbSim.babundance;")
println("...babundance complete...")
for (run_id,) in execute(dbSim, "SELECT run_id FROM runs ORDER BY run_id")
    println(run_id)
    virusesDF = DataFrame(execute(
        dbSim,
        "SELECT t, viral_abundance FROM summary WHERE run_id = $(run_id) AND viral_abundance = 0 ORDER by t DESC"
    ))
    if size(virusesDF)[1] != 0
        (tExt,) = execute(
            dbSim,
            "SELECT t_extinction FROM vextinctions 
            WHERE run_id = $(run_id) ORDER BY t_extinction DESC LIMIT 1")
        execute(dbOutput, "INSERT INTO vextinctions VALUES ($(tExt[1]),$(run_id))")
    end
end
execute(dbOutput, "COMMIT")

println("Complete!")