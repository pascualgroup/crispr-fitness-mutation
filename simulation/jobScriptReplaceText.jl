#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
include("src/setup.jl")

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
db = SQLite.DB(joinpath(SCRIPT_PATH, "sweep_db.sqlite"))

function replaceText()
    for (jobID,) in execute(db, "SELECT job_id FROM jobs ORDER BY job_id")
        jobPath = joinpath(SCRIPT_PATH, "jobs","$(jobID)","job.sbatch")
        run(`sed -i '/#SBATCH --time=1-12:00:00/c\#SBATCH --time=9:00:00' $(jobPath)`)
        println("job script $(jobID) complete") 
    end
end

replaceText()

#db = SQLite.DB(joinpath("sweep_db_gathered.sqlite"))
#runIDs = [runID for (runID,) in execute(db, "SELECT DISTINCT run_id FROM summary ORDER BY run_id")]
