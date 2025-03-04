#!/usr/bin/env julia

"""
This script gathers the output generated by the sweep in `generate-sweep.jl`
into a single SQLite database.
This is probably more than you want, but it shows generally how to consolidate
things.
"""

println("(Julia compilation delay...)")

include(joinpath("model","setup.jl"))

using SQLite

## Define Paths ##
analysisType = "$(ARGS[1])"
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
dbSimPath = joinpath(SCRIPT_PATH,"..","simulation","sweep_db_gathered.sqlite") # cluster
analysisDir = joinpath(SCRIPT_PATH,"analyses",analysisType) # cluster
##

if length(ARGS[2:end]) == 0
    error("please include more info as arguments")
end

if analysisType == "tripartite-networks"
    if !in(ARGS[2],["plots","plot","data"])
        error("`please indicate the second argument as `data` or `plots`")
    end
    println("Please enter the run ids to be analyzed; press Ctrl-D when done.")
    ids_entered = readlines()
    run_ids = map(x->parse(Int64,x), ids_entered)
    if ARGS[2] == "plots" || ARGS[2] == "plot"
        if analysisType == "matches"
            error("plots for `matches.jl` does not exist")
        end
        ROOT_RUN_SCRIPT = joinpath(SCRIPT_PATH,"analyses","networks","make-$(analysisType)-plots.py")
        # println(ROOT_RUN_SCRIPT)
    elseif ARGS[2] == "data"
        ROOT_RUN_SCRIPT = joinpath(SCRIPT_PATH,"analyses","networks","$(analysisType).jl")
    else
        error("invalid input")
    end
elseif analysisType == "walls-shannon"
    if !in(ARGS[2],["plots","plot","data"])
        error("`please indicate the second argument as `data` or `plots`")
    end
    if length(ARGS[3:end]) < 3
        error("at least 1 of 3 wall threshold values are missing")
    end
    if analysisType == "walls-shannon" && parse(Float64,ARGS[2]) <= parse(Float64,ARGS[3])
        error("upper percent threshold must be greater than lower percent threshold")
    end

    if analysisType == "walls-shannon" && parse(Float64,ARGS[2]) > 100
        error("upper percent cannot be greater than 100%")
    end
    println("Please enter the run ids to be analyzed; press Ctrl-D when done.")
    ids_entered = readlines()
    run_ids = map(x->parse(Int64,x), ids_entered)
    if ARGS[2] == "plots" || ARGS[2] == "plot"
        ROOT_RUN_SCRIPT = joinpath(analysisDir,"make-$(analysisType)-plots.py")
    elseif ARGS[2] == "data"
        ROOT_RUN_SCRIPT = joinpath(analysisDir,"$(analysisType).jl")
    else
        error("invalid input")
    end
else
    if !in(ARGS[2],["plots","plot","data"])
        error("`please indicate the second argument as `data` or `plots`")
    end
    println("Please enter the run ids to be analyzed; press Ctrl-D when done.")
    ids_entered = readlines()
    run_ids = map(x->parse(Int64,x), ids_entered)
    if ARGS[2] == "plots" || ARGS[2] == "plot"
        if analysisType == "matches"
            error("plots for `matches.jl` does not exist")
        end
        ROOT_RUN_SCRIPT = joinpath(analysisDir,"make-$(analysisType)-plots.py")
        # println(ROOT_RUN_SCRIPT)
    elseif ARGS[2] == "data"
        ROOT_RUN_SCRIPT = joinpath(analysisDir,"$(analysisType).jl")
    else
        error("invalid input")
    end
end



ROOT_RUNMANY_SCRIPT = joinpath(SCRIPT_PATH,"model", "runmany.jl")
cd(SCRIPT_PATH)

# Number of SLURM jobs to generate
const N_JOBS_MAX = 1
const N_CORES_PER_JOB_MAX = 1 # Half a node (14) is easier to get scheduled than a whole one
const mem_per_cpu = 8000 # in MB 1000MB = 1 GB

function main()
    # Root run directory
    if isfile(joinpath("isolates","$(analysisType)jobs.sqlite"))
        rm(joinpath("isolates","$(analysisType)jobs.sqlite"),force=true)
        println("deleted `/isolates/$(analysisType)jobs.sqlite`.")
    end

    if ispath(joinpath("isolates","$(analysisType)-jobs"))
        rm(joinpath("isolates","$(analysisType)-jobs"),force=true,recursive=true)
        println("deleted `/isolates/$(analysisType)-jobs`.")
    end
    if !ispath(joinpath("isolates"))
        mkpath(joinpath("isolates"))
    end
    dbSim = SQLite.DB(dbSimPath)
    # Create little database that corresponds analysis runs to jobIDs for troubleshooting
    dbTempJobs = SQLite.DB(joinpath("isolates","$(analysisType)jobs.sqlite"))
    execute(dbTempJobs, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(dbTempJobs, "CREATE TABLE job_runs (job_id INTEGER, run_id INTEGER, run_dir TEXT)")

    numSubmits = generate_analysis_runs(dbSim)
    generate_analysis_jobs(dbSim,dbTempJobs,numSubmits)
end

function generate_analysis_runs(dbSim::DB) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    run_count = 0
    println("Processing analysis script for each run")
    for run_id in run_ids
        (cr,) = execute(dbSim, "SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
        run_dir = joinpath("isolates", "runID$(run_id)-c$(cr[1])-r$(cr[2])")
        if ARGS[2] == "data"
            if isfile(joinpath(run_dir,"$(analysisType)_output.sqlite"))
                error("Please move, delete, or rename `/isolates/runID$(run_id)-c$(cr[1])-r$(cr[2])/$(analysisType)_output.sqlite`.")
            end
        end
        if !ispath(run_dir)
            mkpath(run_dir)
        end

        # Generate shell script to perform a single run
        run_script = joinpath(run_dir, "run-$(analysisType).sh")
        if length(ARGS[3:end]) > 0
            if ARGS[2] == "data"
                argString = map(x->string("$(x) "), ARGS[3:end]) # the space after $(x) is important
                open(run_script, "w") do f
                    print(f, """
                    #!/bin/sh
                    cd `dirname \$0`
                    julia $(ROOT_RUN_SCRIPT) $(run_id) $(argString...) &> run-output-$(analysisType).txt
                    """)
                end
            end
            if in(ARGS[2], ["plot","plots"])
                argString = map(x->string("$(x) "), ARGS[3:end]) # the space after $(x) is important
                open(run_script, "w") do f
                    print(f, """
                    #!/bin/sh
                    cd `dirname \$0`
                    module load python/anaconda-2021.05
                    python $(ROOT_RUN_SCRIPT) $(run_id) $(argString...) &> run-output-$(analysisType).txt
                    """)
                end
            end
        else
            if ARGS[2] == "data"
                open(run_script, "w") do f
                    print(f, """
                    #!/bin/sh
                    cd `dirname \$0`
                    julia $(ROOT_RUN_SCRIPT) $(run_id) &> run-output-$(analysisType).txt
                    """)
                end
            end
            if in(ARGS[2], ["plot","plots"])
                open(run_script, "w") do f
                    print(f, """
                    #!/bin/sh
                    cd `dirname \$0`
                    module load python/anaconda-2021.05
                    python $(ROOT_RUN_SCRIPT) $(run_id) &> run-output-$(analysisType).txt
                    """)
                end
            end
        end
        run(`chmod +x $(run_script)`) # Make run script executable
        run_count += 1
    end
    return numSubmits = Int64(ceil(run_count/(N_JOBS_MAX*N_CORES_PER_JOB_MAX)))
end

function generate_analysis_jobs(dbSim::DB,dbTempJobs::DB,numSubmits::Int64)
    println("Assigning analysis runs to jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    n_cores_count = 0

    execute(dbTempJobs, "BEGIN TRANSACTION")

    for run_id in run_ids
        (cr,) = execute(dbSim, "SELECT combo_id,replicate FROM runs WHERE run_id = $(run_id)")
        run_dir = joinpath("isolates", "runID$(run_id)-c$(cr[1])-r$(cr[2])")
        execute(dbTempJobs, "INSERT INTO job_runs VALUES (?,?,?)", (job_id, run_id, run_dir))
        # Mod-increment job ID
        job_id = mod(job_id,N_JOBS_MAX*numSubmits) + 1
    end

    submitScripts = IOStream[]
    for script in 1:numSubmits
        push!(submitScripts,open("$(script)_$(analysisType)-submit-jobs.sh", "w"))
        println(submitScripts[script], """
        #!/bin/sh
        cd `dirname \$0`
        """)
    end

    for (job_id,) in execute(dbTempJobs, "SELECT DISTINCT job_id FROM job_runs ORDER BY job_id")

        job_dir = joinpath("isolates","$(analysisType)-jobs", "$(job_id)")
        @assert !ispath(job_dir)
        mkpath(job_dir)

        # Get all run directories for this job
        run_dirs = [run_dir for (run_dir,) in execute(dbTempJobs,
            """
            SELECT run_dir FROM job_runs
            WHERE job_id = ?
            """,
            (job_id,)
        )]

        n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX)

        if n_cores > n_cores_count
            n_cores_count = n_cores
        end

        # Write out list of runs
        open(joinpath(job_dir, "runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(SCRIPT_PATH, run_dir, "run-$(analysisType).sh")
                println(f, run_script)
            end
        end

        # Create job sbatch file
        job_sbatch = joinpath(job_dir, "job.sbatch")
        open(job_sbatch, "w") do f
            print(f, """
            #!/bin/sh
            #SBATCH --account=pi-pascualmm
            #SBATCH --partition=broadwl
            #SBATCH --job-name=$(job_id)I$(analysisType)
            #SBATCH --tasks=1
            #SBATCH --cpus-per-task=$(n_cores)
            #SBATCH --mem-per-cpu=$(mem_per_cpu)m
            #SBATCH --time=1-12:00:00
            #SBATCH --chdir=$(joinpath(SCRIPT_PATH, job_dir))
            #SBATCH --output=job-output.txt
            #SBATCH --mail-user=armun@uchicago.edu
            # Uncomment this to use the Midway-provided Julia:
            module load julia
            julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) runs.txt
            """)
        end
        run(`chmod +x $(job_sbatch)`) # Make run script executable (for local testing)

        execute(dbTempJobs, "INSERT INTO jobs VALUES (?,?)", (job_id, job_dir))

        submitScript = mod(job_id,numSubmits) + 1
        println(submitScripts[submitScript], "sbatch $(job_sbatch)")
    end
    execute(dbTempJobs, "COMMIT")
    map(close,submitScripts)

    #run(`chmod +x submit_analysis_jobs.sh`) # Make submit script executable
    @info "
    Sweep will be submitted via $(numSubmits) `$(analysisType)-submit-jobs.sh` script(s).
    Each `analysis_submit_jobs.sh` script submits $(N_JOBS_MAX) job(s).
    Each job will use $(n_cores_count) cpus (cores) at most, where each cpu will use $(mem_per_cpu/1000)GB.
    Each job therefore will use at most $(n_cores_count*mem_per_cpu/1000)GB of memory in total.
    "
end


main()
