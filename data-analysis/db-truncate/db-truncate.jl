#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]
##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah","crispr-sweep-7-2-2022/isolates/runID3297-c66-r47/runID3297-c66-r47.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/escape-count.sqlite") # local
dbOutputPath = joinpath("db-truncate_output.sqlite") # cluster
rm(dbOutputPath, force=true)
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)

execute(dbOutput, "CREATE TABLE bgrowthrates (bstrain_id INTEGER,  growth_rate REAL)")
execute(dbOutput, "CREATE TABLE summary (t INTEGER, viral_abundance INTEGER,  microbial_abundance INTEGER)")
execute(dbOutput, "CREATE TABLE babundance (t INTEGER, bstrain_id INTEGER,  abundance INTEGER)")
execute(dbOutput, "CREATE TABLE vextinctions (t_extinction REAL)")


execute(dbOutput, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbOutput, "BEGIN TRANSACTION")
execute(dbOutput, "INSERT INTO bgrowthrates (bstrain_id, growth_rate) SELECT bstrain_id, growth_rate FROM dbSim.bgrowthrates WHERE run_id = $(run_id);")
execute(dbOutput, "INSERT INTO summary (t, viral_abundance, microbial_abundance) SELECT t, viral_abundance, microbial_abundance 
                    FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbOutput, "INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
virusesDF = DataFrame(execute(
    dbSim,
    "SELECT t, viral_abundance FROM summary WHERE run_id = $(run_id) AND viral_abundance = 0 ORDER by t DESC"
))
if size(virusesDF)[1] != 0
    (tExt,) = execute(
        dbSim,
        "SELECT t_extinction FROM vextinctions 
        WHERE run_id = $(run_id) ORDER BY t_extinction DESC LIMIT 1")
    execute(dbOutput, "INSERT INTO vextinctions VALUES ($(tExt[1]))")
end
execute(dbOutput, "COMMIT")

println("Complete!")