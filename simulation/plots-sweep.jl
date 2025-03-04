#!/usr/bin/env julia

"""
The purpose of this file is to illustrate how to do a parameter sweep in Julia.
This script loops through parameter combinations, and replicates with different
random seeds, and generates files necessary to perform runs on a local machine
or on a SLURM cluster.
To use this script for an experiment, you should copy this directory to a new
location, modify the parameter sweeps, and modify the relative paths to
`preamble.jl` and `ROOT_PATH` below to be correct, and the run
`./generate-runs.jl`
When the experiment is complete, you can collect all output files into a single
SQLite database via
`./gather-output.jl`
For each run, it creates a directory, `runs/c<combo_id>/r<replicate>`, and adds
entries to a SQLite database of run information, to make it easy to identify
runs and collate output.
It also divides runs into jobs suitable for execution on a single cluster node
or local machine. The runs are specified as lines in the job's `runs.txt`
file, and the job is specified in a `job.sbatch` file, which can be run directly
as a shell script or submitted to a SLURM cluster.
Each job uses the script `varmodel3/runmany.jl` to run a single-node, multi-core
queue of runs, with one run running on each core at any time.
This script also generates a script `submit_jobs.sh`, which submits every job to
SLURM at once.
Runs are divided into at most `N_JOBS_MAX` jobs that make use of at most
`N_CORES_PER_JOB_MAX` for the cluster node's local queue.
This allows you to work within limits set by your cluster administrator.
If you have no limits, you should set `N_JOBS_MAX` to a very large number,
and set `N_CORES_PER_JOB_MAX = 1`, so that the cluster can dynamically
balance runs across cluster nodes as the experiment runs.
To modify configuration settings for SLURM jobs, edit the template string in
the `generate_jobs()` function.
"""

println("(Julia compilation delay...)")

include("src/setup.jl")

using SQLite

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
# ROOT_PATH = joinpath(SCRIPT_PATH)
ROOT_PATH = joinpath("/scratch/cgsb/pascual/armun/crispr/vary-fitness/20_MOI3/simulation")
# dbSimPath = joinpath(SCRIPT_PATH, "..", "simulation", "sweep_db.sqlite")
dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/20_MOI3/sweep_db.sqlite")
ROOT_RUN_SCRIPT = joinpath(ROOT_PATH, "src","make-plots.py")
ROOT_RUNMANY_SCRIPT = joinpath(ROOT_PATH,"src", "runmany.jl")
cd(SCRIPT_PATH)

# Number of SLURM jobs to generate
const N_JOBS_MAX = 500
const N_PROCESSES = 15 # Half a node, easier to get scheduled than a whole one
const N_CORES_PER_JOB_MAX = 15 # Half a node, easier to get scheduled than a whole one
const mem_per_cpu = 8000 # in MB 1000MB = 1 GB



function main()
    # Root run directory
    # if ispath(joinpath(SCRIPT_PATH,"plots"))
    #     error("`/simulation/plots` already exists; please delete first.")
    # end

    # if !ispath(joinpath(SCRIPT_PATH,"runs"))
    #     error("`/simulation/runs` is missing; please run generate-sweep.jl first, then simulate time series.")
    # end

    # if !isfile(joinpath(SCRIPT_PATH,"sweep_db.sqlite"))
    #     error("`sweep_db.sqlite` is missing; please run generate-sweep.jl first, then simulate time series.")
    # end

    # # Root job directory
    # if !ispath(joinpath(SCRIPT_PATH,"jobs"))
    #     @info "Note that `/simulation/jobs` is missing."
    # end

    dbSim = SQLite.DB(dbSimPath)
    db = SQLite.DB(joinpath("plots_db.sqlite"))
    execute(db, "DROP TABLE IF EXISTS plot_jobs")
    execute(db, "CREATE TABLE plot_jobs (job_id INTEGER, job_dir TEXT)")
    execute(db, "DROP TABLE IF EXISTS plot_job_runs")
    execute(db, "CREATE TABLE plot_job_runs (job_id INTEGER, run_id INTEGER, run_dir TEXT)")

    numSubmits = generate_plot_runs(dbSim,db)
    generate_plot_jobs(dbSim,db, numSubmits)
end

function generate_plot_runs(dbSim::DB,db::DB) # This function generates the directories
    # for the individual parameter sets and corresponding replicates. It also
    # generates shell scripts for each run and corresponding parameter file.

    # Loop through parameter combinations and replicates, generating a run directory
    # `plots/c<combo_id>/r<replicate>` for each one.
    run_count = 0
    println("Processing plot script for each run")
    for (run_id, combo_id, replicate) in execute(dbSim, "SELECT run_id,combo_id,replicate FROM runs ORDER BY combo_id,replicate")
        #println("Processing plot script for combination $(combo_id)/replicate $(replicate)"
        #) # local
        # if run_id in [1, 101, 201, 7, 107, 207, 8, 108, 208]
        #     println("...skipping run_id $(run_id)...")
        #     continue
        # end
        # run_dir = joinpath(SCRIPT_PATH,"runs", "c$(combo_id)", "r$(replicate)")
        # @assert ispath(run_dir)
        run_dir = joinpath(SCRIPT_PATH, "plots", "runs", "c$(combo_id)", "r$(replicate)")

        plot_dir = joinpath(SCRIPT_PATH,"plots", "runs","c$(combo_id)", "r$(replicate)")
        # @assert !ispath(plot_dir)
        mkpath(plot_dir)

        argString = map(x->string("$(x) "), ARGS) # the space after $(x) is important

        # Generate shell script to perform a single run
        run_script = joinpath(run_dir, "runplotmaker.sh")
        open(run_script, "w") do f
            print(f, """
            #!/bin/sh
            cd `dirname \$0`
            module load python/intel/3.8.6
            python $(ROOT_RUN_SCRIPT) $(run_id) $(argString...) &> plot_output.txt
            """)
        end
        run(`chmod +x $(run_script)`)
        run_count += 1
    end
    return numSubmits = Int64(ceil(run_count/(N_JOBS_MAX*N_PROCESSES)))
end

function generate_plot_jobs(dbSim::DB,db::DB,numSubmits::Int64)
    println("Assigning plot runs to plot jobs...")

    # Assign runs to jobs (round-robin)
    job_id = 1
    n_cores = 0

    execute(db, "BEGIN TRANSACTION")

    for (run_id, run_dir) in execute(dbSim, "SELECT run_id, run_dir FROM runs ORDER BY combo_id, replicate")
        execute(db, "INSERT INTO plot_job_runs VALUES (?,?,?)", (job_id, run_id, joinpath("plots",run_dir)))
        # Mod-increment job ID
        job_id = mod(job_id,N_JOBS_MAX*numSubmits) + 1
    end

    submitScripts = IOStream[]
    for script in 1:numSubmits
        push!(submitScripts,open("$(script)_plot-submit-jobs.sh", "w"))
        println(submitScripts[script], """
        #!/bin/sh
        cd `dirname \$0`
        """)
    end

    for (job_id,) in execute(db, "SELECT DISTINCT job_id FROM plot_job_runs ORDER BY job_id")

        job_dir = joinpath("plots","jobs", "$(job_id)")
        mkpath(job_dir)
        @assert ispath(job_dir)

        # Get all run directories for this job
        run_dirs = [run_dir for (run_dir,) in execute(db,
            """
            SELECT run_dir FROM plot_job_runs
            WHERE job_id = ?
            """,
            (job_id,)
        )]


        if N_CORES_PER_JOB_MAX == nothing
            n_cores = min(length(run_dirs), N_PROCESSES)
        else
            n_cores = N_CORES_PER_JOB_MAX
        end

        # Write out list of runs
        open(joinpath(job_dir, "plot_runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(ROOT_PATH, run_dir, "runplotmaker.sh")
                println(f, run_script)
            end
        end

        # Create job sbatch file
        job_sbatch = joinpath(job_dir, "plotjob.sbatch")
        open(job_sbatch, "w") do f
            print(f, """
            #!/bin/sh
            #SBATCH --nodes=1
            #SBATCH --ntasks-per-node=$(n_cores)
            #SBATCH --cpus-per-task=1
            #SBATCH --job-name=$(job_id)plots
            #SBATCH --mem-per-cpu=$(mem_per_cpu)m
            #SBATCH --time=01:55:00
            #SBATCH --chdir=$(joinpath(ROOT_PATH, job_dir))
            #SBATCH --output=output.txt
            #SBATCH --mail-type=END,FAIL
            #SBATCH --mail-user=al8784@nyu.edu
            module purge
            module load julia/1.6.1
            ~/julia/my-julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) plot_runs.txt
            """) # runs.txt is for parallel processing
        end
        run(`chmod +x $(job_sbatch)`) # Make run script executable (for local testing)

        execute(db, "INSERT INTO plot_jobs VALUES (?,?)", (job_id, job_dir))

        submitScript = mod(job_id,numSubmits) + 1
        println(submitScripts[submitScript], "sbatch $(job_sbatch)")
    end
    execute(db, "COMMIT")
    map(close,submitScripts)

    #run(`chmod +x submit_plot_jobs.sh`) # Make submit script executable
    @info "
    Sweep will be submitted via $(numSubmits) `plot_submit_jobs.sh` script(s).
    Each `plot_submit_jobs.sh` script submits $(N_JOBS_MAX) jobs.
    Each job will use $(n_cores) cpus (cores) at most, where each cpu will use $(mem_per_cpu/1000)GB.
    Each job therefore will use at most $(n_cores*mem_per_cpu/1000)GB of memory in total.
    "
end


main()
