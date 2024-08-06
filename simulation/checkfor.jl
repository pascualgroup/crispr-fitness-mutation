#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
include("src/setup.jl")

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
db = SQLite.DB(joinpath(SCRIPT_PATH, "sweep_db.sqlite"))

function checkRuns()
    for (run_dir,) in execute(db, "SELECT run_dir FROM runs ORDER BY combo_id, replicate")
        if ispath(joinpath(SCRIPT_PATH, "$(run_dir)","output.txt"))
            check = let x = 0
                for line in eachline(open(joinpath(SCRIPT_PATH, "$(run_dir)","output.txt")))
                    if contains(line,"initial output")
                        x += 1
                    end
                    if contains(line, "2000.0")
                        x += 1
                    end
                end
                x
            end
            if check != 2
                println("error: $(run_dir)")
            end
        else
            println("does not exist: $(run_dir)")
        end
        check = 0
    end
end

function checkJobs()
    check = 0
    for (job_dir,) in execute(db, "SELECT job_dir FROM jobs ORDER BY job_id")
        if ispath(joinpath(SCRIPT_PATH, "$(job_dir)","runs.txt"))
            for line in eachline(open(joinpath(SCRIPT_PATH, "$(job_dir)","runs.txt")))
                if contains(line,"/project2/pascualmm/armun/crispr/vary-fitness/12_MOI3/simulation/runs/c39/r103/run.sh")
                    check += 1
                end
            end
            if check >= 1
                println("found: $(job_dir)")
            end
        else
            println("does not exist: $(job_dir)")
        end
        check = 0
    end
end

checkRuns()
checkJobs()

#db = SQLite.DB(joinpath("sweep_db_gathered.sqlite"))
#runIDs = [runID for (runID,) in execute(db, "SELECT DISTINCT run_id FROM summary ORDER BY run_id")]
