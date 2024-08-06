
### PARAMETERS ###

@with_kw mutable struct Params
    "Simulation end time"
    t_final::Union{Float64,Nothing}

    "Time between output events"
    t_output::Union{Float64,Nothing}

    "Seed for random number generator"
    rng_seed::Union{UInt64,Nothing}

    "Seed for random number generator for community initialization"
    rngDiv_seed::Union{UInt64,Nothing}

    "Enable output?"
    enable_output::Union{Bool,Nothing}

    "Truncate output based on viral extinction?"
    truncate_output::Union{Bool,Nothing}

    "Auxiliary parameter for miscellaneous use"
    auxiliary::Any

    "Choose type of microbe initial community condition"
    init_bcomm_function::Union{UInt64,Nothing}

    "Choose type of viral initial community condition"
    init_vcomm_function::Union{UInt64,Nothing}

    "Directional growth mutations?"
    directional::Union{Bool,Nothing}

    "The is the genotype to phenotype map for microbial demographic trait"
    evofunction::Union{UInt64,Nothing}

    "This is an arbitrary scale factor used for the genotype to phenotype map: evofunction"
    evofunctionScale::Union{Float64,Nothing}

    "This sets the initial allele value for the host when a varied community is not initialized"
    initial_locus_allele::Union{Bool,Float64,Nothing}

    "This sets the middle allele value"
    center_allele::Union{Float64,Nothing}

    "The amount allele value changes upon directional mutation [epsilon]"
    locus_resolution::Union{Float64,Nothing}

    "This sets the upper and lower bounds of allele values"
    max_allele::Union{Float64,Nothing}

    "This sets maximum intrinsic fitness allowed among microbes"
    max_fitness::Union{Float64,Nothing}

    "Number of initial bacterial strains"
    n_bstrains::Union{UInt64,Nothing}

    "Number of initial hosts per bacterial strain"
    n_hosts_per_bstrain::Union{UInt64,Nothing}

    "Number of initial virus strains"
    n_vstrains::Union{UInt64,Nothing}

    "Number of initial particles per bacterial strain"
    n_particles_per_vstrain::Union{UInt64,Nothing}

    "Number of initial protospacers per virus strain"
    n_protospacers::Union{UInt64,Nothing}

    "Maximum number of spacers in a bacterial strain"
    n_spacers_max::Union{UInt64,Nothing}

    "CRIPSR failure probability [p]"
    crispr_failure_prob::Union{Float64,Nothing}

    "New spacer acquisition probability [q]"
    spacer_acquisition_prob::Union{Float64,Nothing}

    "Microbial mutation probability [rho]"
    microbe_mutation_prob_per_replication::Union{Float64,Nothing}

    "Microbial mutation probability [rho]"
    microbe_mutation_prob_per_spacer::Union{Float64,Nothing}

    "Carrying capacity (1/mL) = [K]"
    microbe_carrying_capacity::Union{Float64,Nothing}

    "Burst size [beta]"
    viral_burst_size::Union{UInt64,Nothing}

    "Adsorption rate [phi]"
    adsorption_rate::Union{Float64,Nothing}

    "Viral decay rate [m]"
    viral_decay_rate::Union{Float64,Nothing}

    "Mutation rate [mu]"
    viral_mutation_rate::Union{Float64,Nothing}

    "Constant death rate (not in Childs model) [d]"
    microbe_death_rate::Union{Float64,Nothing}

    "Constant immigration rate (not in Childs model) [eta]"
    microbe_immigration_rate::Union{Float64,Nothing}
end

@with_kw mutable struct ParamSweep
    "Simulation end time"
    t_final::Union{Float64,Array{Float64,1},Nothing}

    "Time between output events"
    t_output::Union{Float64,Array{Float64,1},Nothing}

    "Seed for random number generator"
    rng_seed::Union{UInt64,Array{UInt64,1},Nothing}

    "Seed for random number generator for community initialization"
    rngDiv_seed::Union{UInt64,Array{UInt64,1},Nothing}

    "Enable output?"
    enable_output::Union{Bool,Nothing}

    "Truncate output based on viral extinction?"
    truncate_output::Union{Bool,Nothing}

    "Auxiliary parameter for miscellaneous use"
    auxiliary::Any

    "Choose type of microbe initial community condition"
    init_bcomm_function::Union{UInt64,Array{UInt64,1},Nothing}

    "Choose type of viral initial community condition"
    init_vcomm_function::Union{UInt64,Array{UInt64,1},Nothing}

    "Directional growth mutations?"
    directional::Union{Bool,Nothing}

    "The is the genotype to phenotype map for microbial demographic trait"
    evofunction::Union{UInt64,Nothing}

    "This is an arbitrary scale factor used for the genotype to phenotype map: evofunction"
    evofunctionScale::Union{Float64,Array{Float64,1},Nothing}

    "This sets the initial allele value"
    initial_locus_allele::Union{Bool,Float64,Array{Union{Bool,Float64},1},Nothing}

    "This sets the middle allele value"
    center_allele::Union{Float64,Array{Float64,1},Nothing}

    "The amount allele value changes upon directional mutation [epsilon]"
    locus_resolution::Union{Float64,Array{Float64,1},Nothing}

    "This sets the upper and lower bounds of allele values"
    max_allele::Union{Float64,Array{Float64,1},Nothing}

    "This sets the maximum instrinsic fitness of microbe"
    max_fitness::Union{Float64,Array{Float64,1},Nothing}

    "Number of initial bacterial strains"
    n_bstrains::Union{UInt64,Array{UInt64,1},Nothing}

    "Number of initial hosts per bacterial strain"
    n_hosts_per_bstrain::Union{UInt64,Array{UInt64,1},Nothing}

    "Number of initial virus strains"
    n_vstrains::Union{UInt64,Array{UInt64,1},Nothing}

    "Number of initial particles per bacterial strain"
    n_particles_per_vstrain::Union{UInt64,Array{UInt64,1},Nothing}

    "Number of initial protospacers per virus strain"
    n_protospacers::Union{UInt64,Array{UInt64,1},Nothing}

    "Maximum number of spacers in a bacterial strain"
    n_spacers_max::Union{UInt64,Array{UInt64,1},Nothing}

    "CRIPSR failure probability [p]"
    crispr_failure_prob::Union{Float64,Array{Float64,1},Nothing}

    "New spacer acquisition probability [q]"
    spacer_acquisition_prob::Union{Float64,Array{Float64,1},Nothing}

    "Microbial mutation probability [rho]"
    microbe_mutation_prob_per_replication::Union{Float64,Array{Float64,1},Nothing}

    "Microbial mutation probability [rho]"
    microbe_mutation_prob_per_spacer::Union{Float64,Array{Float64,1},Nothing}

    "Carrying capacity (1/mL) = [K]"
    microbe_carrying_capacity::Union{Float64,Array{Float64,1},Nothing}

    "Burst size [beta]"
    viral_burst_size::Union{UInt64,Array{UInt64,1},Nothing}

    "Adsorption rate [phi]"
    adsorption_rate::Union{Float64,Array{Float64,1},Nothing}

    "Viral decay rate [m]"
    viral_decay_rate::Union{Float64,Array{Float64,1},Nothing}

    "Mutation rate [mu]"
    viral_mutation_rate::Union{Float64,Array{Float64,1},Nothing}

    "Constant death rate (not in Childs model) [d]"
    microbe_death_rate::Union{Float64,Array{Float64,1},Nothing}

    "Constant immigration rate (not in Childs model) [eta]"
    microbe_immigration_rate::Union{Float64,Array{Float64,1},Nothing}
end


function init_params(d_symb::Dict{Symbol,Any})
    Params(;
        t_final = d_symb[:t_final][1],

        t_output = d_symb[:t_output][1],

        rng_seed = d_symb[:rng_seed][1],

        rngDiv_seed=d_symb[:rngDiv_seed][1],

        enable_output = d_symb[:enable_output][1],

        truncate_output=d_symb[:truncate_output][1],

        auxiliary = d_symb[:auxiliary][1],

        init_bcomm_function = d_symb[:init_bcomm_function][1],

        init_vcomm_function=d_symb[:init_vcomm_function][1],

        directional = d_symb[:directional][1],

        evofunction = d_symb[:evofunction][1],

        evofunctionScale = d_symb[:evofunctionScale][1],

        initial_locus_allele = d_symb[:initial_locus_allele][1],

        center_allele = d_symb[:center_allele][1],

        locus_resolution = d_symb[:locus_resolution][1],

        max_allele = d_symb[:max_allele][1],

        max_fitness = d_symb[:max_fitness][1],

        n_bstrains = d_symb[:n_bstrains][1],

        n_hosts_per_bstrain = d_symb[:n_hosts_per_bstrain][1],

        n_vstrains = d_symb[:n_vstrains][1],

        n_particles_per_vstrain = d_symb[:n_particles_per_vstrain][1],

        n_protospacers = d_symb[:n_protospacers][1],

        n_spacers_max = d_symb[:n_spacers_max][1],

        crispr_failure_prob = d_symb[:crispr_failure_prob][1],

        spacer_acquisition_prob = d_symb[:spacer_acquisition_prob][1],

        microbe_mutation_prob_per_replication=d_symb[:microbe_mutation_prob_per_replication][1],

        microbe_mutation_prob_per_spacer=d_symb[:microbe_mutation_prob_per_spacer][1],

        microbe_carrying_capacity = d_symb[:microbe_carrying_capacity][1],

        viral_burst_size = d_symb[:viral_burst_size][1],

        adsorption_rate = d_symb[:adsorption_rate][1],

        viral_decay_rate = d_symb[:viral_decay_rate][1],

        viral_mutation_rate = d_symb[:viral_mutation_rate][1],

        microbe_death_rate = d_symb[:microbe_death_rate][1],

        microbe_immigration_rate = d_symb[:microbe_immigration_rate][1],
    )
end

### SIMULATION STATE ###

mutable struct Strains
    next_id::UInt64
    ids::Vector{UInt64}

    abundance::Vector{UInt64}
    total_abundance::UInt64

    spacers::Vector{Vector{UInt64}}
    growthalleles::Vector{Float64}
    growthrates::Vector{Float64}
end


# function make_bstrains(fxnNum, rng, rngDiv, n_strains, n_hosts_per_strain, 
#     vstrainIDs, vpspacers,
#     microbe_carrying_capacity,
#     initial_locus_allele, evofunction,
#     maxallele, centerallele, scale, maxfitness, death)
#     if initialConditionsDB == Nothing
#         initialalleles = repeat([initial_locus_allele], n_strains)
#         initialgrowthrates = map(x -> fitness(x, centerallele, maxallele,
#                 scale, maxfitness, death, evofunction), initialalleles)
#         ids = Vector(1:n_strains)
#         abundance = repeat([n_hosts_per_strain], n_strains)
#         total_abundance = n_strains * n_hosts_per_strain
#         Strains(
#             n_strains + 1, # without initial conditions we just have one strain, and the next_strain_id is just 2
#             ids,
#             abundance,
#             total_abundance,
#             repeat([[]], n_strains),
#             initialalleles,
#             initialgrowthrates
#         )
#     else
#         ids = [strain_id for (strain_id,) in
#                execute(initialConditionsDB, "SELECT bstrain_id FROM babundance ORDER BY bstrain_id")]
#         abundance = [abund for (abund,) in
#                      execute(initialConditionsDB, "SELECT abundance FROM babundance ORDER BY bstrain_id")]
#         if random_initial_alleles
#             locusalleles = [centerallele, 
#                             rand(rng, Uniform(centerallele - maxallele, maxallele), 
#                                 length(ids)-1)...]
#         else
#             locusalleles = [allele for (allele,) in
#                             execute(initialConditionsDB, "SELECT locus_allele FROM blocusalleles ORDER BY bstrain_id")]
#         end
#         growthrates = map(x -> fitness(x, centerallele, maxallele,
#                 scale, maxfitness, death, evofunction), locusalleles)
#         total_abundance = sum(abundance)
#         spacers = [
#             Vector([spacer_id for (spacer_id,) in
#                     execute(
#                 initialConditionsDB,
#                 "SELECT spacer_id FROM
#                 bspacers WHERE bstrain_id = $(id)
#                 ORDER BY spacer_id"
#             )])
#             for id in ids
#         ]
#         spacers[spacers.==[[0]]] = repeat([[]], length(spacers[spacers.==[[0]]]))
#         Strains(
#             maximum(ids) + 1,
#             ids,
#             abundance,
#             total_abundance,
#             spacers,
#             locusalleles,
#             growthrates
#         )
#     end
# end

# function make_vstrains(fxnNum, rng, rngDiv, n_strains, n_particles_per_strain, n_pspacers_init)
#     if initialConditionsDB == Nothing
#         next_id = UInt64(n_strains + 1)
#         # println("next id is $(next_id)")
#         # println("type is $(typeof(next_id))")
#         ids = Vector(1:n_strains)
#         abundance = repeat([n_particles_per_strain], n_strains)
#         total_abundance = n_strains * n_particles_per_strain
#         pspacers = [
#             Vector(1:n_pspacers_init) .+ repeat([n_pspacers_init * (i - 1)], n_pspacers_init)
#             for i = 1:n_strains
#         ]

#         Strains(
#             next_id,
#             ids,
#             abundance,
#             total_abundance,
#             pspacers,
#             repeat([], n_strains),
#             repeat([], n_strains)
#         )
#     else
#         ids = [strain_id for (strain_id,) in
#                execute(initialConditionsDB, "SELECT vstrain_id FROM vabundance ORDER BY vstrain_id")]
#         abundance = [abund for (abund,) in
#                      execute(initialConditionsDB, "SELECT abundance FROM vabundance ORDER BY vstrain_id")]
#         total_abundance = sum(abundance)
#         pspacers = [
#             Vector([spacer_id for (spacer_id,) in
#                     execute(
#                 initialConditionsDB,
#                 "SELECT spacer_id FROM
#                 vpspacers WHERE vstrain_id = $(id)
#                 ORDER BY spacer_id"
#             )])
#             for id in ids
#         ]
#         Strains(
#             maximum(ids) + 1,
#             ids,
#             abundance,
#             total_abundance,
#             pspacers,
#             repeat([], length(ids)),
#             repeat([], length(ids))
#         )
#     end
# end


#########################################################
##################THIS SAVES NEW STRAINS#################

mutable struct State
    bstrains::Strains
    vstrains::Strains
    next_pspacer_id::UInt64

    function State(
        auxiliary, rngDiv, n_bstrains, n_hosts_per_bstrain, 
        initial_locus_allele, evofunction,
        maxallele, locusresolution, centerallele, scale, maxfitness, death,
        n_vstrains, n_particles_per_vstrain, n_pspacers_init,
        bFxnNum, vFxnNum
    )
        # if initialConditionsDB == Nothing
        #     next_pspacer_id = 1 + n_vstrains * n_pspacers_init
        # else
        next_pspacer_id, vstrains = setVirusCommunity(bFxnNum, auxiliary, n_vstrains, n_particles_per_vstrain, n_pspacers_init)
        # end
        new(
            setMicrobeCommunity(bFxnNum, auxiliary, rngDiv, n_bstrains, n_hosts_per_bstrain,
                vstrains.spacers, 
                initial_locus_allele,
                evofunction, maxallele, locusresolution, centerallele, scale, maxfitness, death),
            vstrains,
            next_pspacer_id
        )
    end
end

function State(p::Params)
    State(
        p.auxiliary,
        MersenneTwister(p.rngDiv_seed), UInt64(p.n_bstrains), UInt64(p.n_hosts_per_bstrain),
        Float64(p.initial_locus_allele),
        p.evofunction, p.max_allele, p.locus_resolution, p.center_allele, p.evofunctionScale, p.max_fitness, p.microbe_death_rate,
        UInt64(p.n_vstrains), UInt64(p.n_particles_per_vstrain), UInt64(p.n_protospacers), 
        UInt64(p.init_bcomm_function), UInt64(p.init_vcomm_function)
    )
end

#########################################################

mutable struct Simulation
    params::Params
    t::Float64
    state::State
    rng::MersenneTwister

    event_rates::Vector{Float64}
    event_counts::Vector{UInt64}

    db::DB

    #meta_file::IOStream
    #summary_file::IOStream

    function Simulation(p::Params)
        #meta_file = open_csv("meta", "key", "value")

        db = initialize_database()

        # Use random seed if provided, or generate one
        rng_seed = p.rng_seed === nothing ? UInt64(rand(RandomDevice(), UInt32)) : p.rng_seed
        p.rng_seed = rng_seed

        #write_csv(meta_file, "rng_seed", rng_seed)

        execute(db,
            "INSERT INTO meta VALUES (?,?)", ["rng_seed", Int64(rng_seed)]
        )
        execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["start_time", Dates.format(start_time, "yyyy-mm-ddTHH:MM:SS")]
        )

        # Initialize & validate model state
        # println("$([p.n_bstrains, p.n_hosts_per_bstrain, p.initial_locus_allele,
        # p.n_vstrains, p.n_particles_per_vstrain, p.n_protospacers])")
        state = State(p)
        validate(state)

        sim = new(
            p, 0.0, state, MersenneTwister(rng_seed),
            zeros(length(EVENTS)), zeros(length(EVENTS)),
            db
        )
        update_rates!(sim)
        sim
    end
end


### VALIDATION TO ENSURE CODE CORRECTNESS IN THE FACE OF OPTIMIZATIONS ###

function validate(p::Union{Params,ParamSweep})
    @assert p.t_final !== nothing
    @assert p.t_output !== nothing

    @assert p.rngDiv_seed !== nothing
    @assert p.directional !== nothing
    @assert p.init_bcomm_function !== nothing
    @assert p.init_vcomm_function !== nothing
    @assert p.evofunction !== nothing
    @assert p.evofunctionScale !== nothing
    @assert p.initial_locus_allele !== nothing
    @assert p.center_allele !== nothing
    @assert p.locus_resolution !== nothing
    @assert p.max_allele !== nothing
    @assert p.max_fitness !== nothing

    @assert p.n_bstrains !== nothing
    @assert p.n_hosts_per_bstrain !== nothing
    @assert p.n_vstrains !== nothing
    @assert p.n_particles_per_vstrain !== nothing
    @assert p.n_protospacers !== nothing
    @assert p.n_spacers_max !== nothing

    @assert p.crispr_failure_prob !== nothing
    @assert p.spacer_acquisition_prob !== nothing
    @assert p.microbe_mutation_prob_per_replication !== nothing
    @assert p.microbe_mutation_prob_per_spacer !== nothing
    @assert p.microbe_carrying_capacity !== nothing
    @assert p.viral_burst_size !== nothing
    @assert p.adsorption_rate !== nothing
    @assert p.viral_decay_rate !== nothing
    @assert p.viral_mutation_rate !== nothing
    @assert p.microbe_death_rate !== nothing
    @assert p.microbe_immigration_rate !== nothing
end

function validate(s::State)
    validate(s.bstrains)
    validate(s.vstrains)
end

function validate(strains::Strains)
    @assert strains.total_abundance == sum(strains.abundance)
    @assert strains.next_id > maximum(strains.ids)
    @assert length(strains.abundance) == length(strains.ids)
    @assert length(strains.spacers) == length(strains.ids)
    if length(strains.growthrates) > 0
        @assert length(strains.growthrates) == length(strains.ids)
    end
end
