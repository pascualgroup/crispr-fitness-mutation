println("(Julia compilation delay...)")

include("simulation/src/setup.jl")
include("simulation/src/structures.jl")
include("simulation/src/util.jl")

SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

function make_sweep_params_script(SCRIPT_PATH)
    numParams = []
    t_final = Float64(2000.0)
    push!(numParams, length(t_final))
    t_output = Float64(1.0)
    push!(numParams, length(t_output))
    rng_seed = UInt64(1)
    push!(numParams, length(rng_seed))
    rngDiv_seed = UInt64(9876543456)
    push!(numParams, length(rngDiv_seed))
    # auxiliary = [["path/to/1", 2, 5.0], ["path/to/2", 2, 5.0]]
    auxiliary = false
    push!(numParams, length(auxiliary))

    enable_output = true
    truncate_output = true
    directional = false
    evofunction = 2
    evofunctionScale = collect(0:0.5:5)
    push!(numParams, length(evofunctionScale))
    init_bcomm_function = [UInt64(1), UInt64(2)]
    push!(numParams, length(init_bcomm_function))
    init_vcomm_function = [UInt64(1)]
    push!(numParams, length(init_vcomm_function))
    initial_locus_allele= Float64(25)
    push!(numParams, length(initial_locus_allele))
    center_allele = Float64(0)
    push!(numParams, length(center_allele))
    locus_resolution = Float64(.02)
    push!(numParams, length(locus_resolution))
    max_allele = Float64(1)
    push!(numParams, length(max_allele))
    max_fitness = Float64(1)
    push!(numParams, length(max_fitness))

    n_bstrains = UInt64(8)
    push!(numParams, length(n_bstrains))
    n_hosts_per_bstrain = UInt64(100)
    push!(numParams, length(n_hosts_per_bstrain))
    n_vstrains = UInt64(8)
    push!(numParams, length(n_vstrains))
    n_particles_per_vstrain = UInt64(100)
    push!(numParams, length(n_particles_per_vstrain))
    n_protospacers = UInt64(15)
    push!(numParams, length(n_protospacers))
    n_spacers_max = UInt64(10)
    push!(numParams, length(n_spacers_max))
    crispr_failure_prob = Float64(0)
    push!(numParams, length(crispr_failure_prob))
    spacer_acquisition_prob = Float64(1.9e-5)
    push!(numParams, length(spacer_acquisition_prob))
    microbe_mutation_prob_per_replication = [Float64(0), Float64(0.25), Float64(0.5), Float64(0.75), Float64(1)]
    push!(numParams, length(microbe_mutation_prob_per_replication))
    microbe_mutation_prob_per_spacer = [Float64(0), Float64(0.25), Float64(0.5), Float64(0.75), Float64(1)]
    push!(numParams, length(microbe_mutation_prob_per_spacer))
    microbe_carrying_capacity = collect(200000:(400000-200000)/5:400000)
    push!(numParams, length(microbe_carrying_capacity))
    viral_burst_size = UInt64(50)
    push!(numParams, length(viral_burst_size))
    adsorption_rate = Float64(1e-07)
    push!(numParams, length(adsorption_rate))
    viral_decay_rate = Float64(0.1)
    push!(numParams, length(viral_decay_rate))
    viral_mutation_rate = Float64(1.04e-06)
    push!(numParams, length(viral_mutation_rate))
    microbe_death_rate = Float64(0.025)
    push!(numParams, length(microbe_death_rate))
    microbe_immigration_rate = Float64(0)
    push!(numParams, length(microbe_immigration_rate))


    params = ParamSweep(;
        t_final = t_final,
        t_output = t_output,
        rng_seed = rng_seed,
        rngDiv_seed = rngDiv_seed,
        enable_output = enable_output,
        truncate_output= truncate_output,
        auxiliary = auxiliary,
        init_bcomm_function = init_bcomm_function,
        init_vcomm_function = init_vcomm_function,
        directional = directional,
        evofunction =  evofunction,
        evofunctionScale =  evofunctionScale,
        initial_locus_allele = initial_locus_allele,
        center_allele = center_allele,
        locus_resolution = locus_resolution,
        max_allele = max_allele,
        max_fitness = max_fitness,
        n_bstrains = n_bstrains,
        n_hosts_per_bstrain = n_hosts_per_bstrain,
        n_vstrains = n_vstrains,
        n_particles_per_vstrain = n_particles_per_vstrain,
        n_protospacers = n_protospacers,
        n_spacers_max = n_spacers_max,
        crispr_failure_prob = crispr_failure_prob,
        spacer_acquisition_prob = spacer_acquisition_prob,
        microbe_mutation_prob_per_replication = microbe_mutation_prob_per_replication,
        microbe_mutation_prob_per_spacer = microbe_mutation_prob_per_spacer,
        microbe_carrying_capacity = microbe_carrying_capacity,
        viral_burst_size = viral_burst_size,
        adsorption_rate = adsorption_rate,
        viral_decay_rate = viral_decay_rate,
        viral_mutation_rate = viral_mutation_rate,
        microbe_death_rate = microbe_death_rate,
        microbe_immigration_rate = microbe_immigration_rate,
    )
    validate(params)
    params_json = pretty_json(params)
    open(joinpath(SCRIPT_PATH, "sweep-parameters.json"), "w") do f
        println(f, params_json)
    end
    println("This sweep has a maximum of $(prod(numParams)) parameter combinations.
    Multiply this value with the number of replicates intended to evaluate numbers of cores necessary.")
end

function pretty_json(params)
    d = Dict(fn => getfield(params, fn) for fn in fieldnames(typeof(params)))
    io = IOBuffer()
    JSON.print(io, d, 2)
    String(take!(io))
end

make_sweep_params_script(SCRIPT_PATH)
