
### PARAMETERS ###

@with_kw mutable struct Params
    "Simulation end time"
    t_final::Union{Float64,Nothing}

    "Time between output events"
    t_output::Union{Float64,Nothing}

    "Seed for random number generator"
    rng_seed::Union{UInt64,Nothing}

    "Enable output?"
    enable_output::Union{Bool,Nothing}

    "Random initial host locus alleles?"
    random_initial_alleles::Union{Bool,Nothing}

    "Directional growth mutations?"
    directional::Union{Bool,Nothing}

    "The is the genotype to phenotype map for microbial demographic trait"
    evofunction::Union{UInt64,Nothing}

    "This is an arbitrary scale factor used for the genotype to phenotype map: evofunction"
    evofunctionScale::Union{Float64,Array{Float64,1},Nothing}

    "This sets the initial allele value for the host when a varied community is not initialized"
    initial_locus_allele::Union{Float64,Nothing}

    "This sets the middle allele value"
    center_allele::Union{Float64,Nothing}

    "The amount allele value changes upon directional mutation [epsilon]"
    allelic_change::Union{Float64,Nothing}

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
    microbe_mutation_prob::Union{Float64,Nothing}

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


function make_bstrains(rng, n_strains, n_hosts_per_strain, random_initial_alleles, 
    initial_locus_allele, evofunction,
    maxallele, centerallele, scale, maxfitness, death,
    initialConditionsDB::Union{DB,Nothing})
    if initialConditionsDB == Nothing
        initialalleles = repeat([initial_locus_allele], n_strains)
        initialgrowthrates = map(x -> fitness(x, centerallele, maxallele,
                scale, maxfitness, death, evofunction), initialalleles)
        ids = Vector(1:n_strains)
        abundance = repeat([n_hosts_per_strain], n_strains)
        total_abundance = n_strains * n_hosts_per_strain
        Strains(
            n_strains + 1, # without initial conditions we just have one strain, and the next_strain_id is just 2
            ids,
            abundance,
            total_abundance,
            repeat([[]], n_strains),
            initialalleles,
            initialgrowthrates
        )
    else
        ids = [strain_id for (strain_id,) in
               execute(initialConditionsDB, "SELECT bstrain_id FROM babundance ORDER BY bstrain_id")]
        abundance = [abund for (abund,) in
                     execute(initialConditionsDB, "SELECT abundance FROM babundance ORDER BY bstrain_id")]
        if random_initial_alleles
            locusalleles = [centerallele, 
                            rand(rng, Uniform(centerallele - maxallele, maxallele), 
                                length(ids)-1)...]
        else
            locusalleles = [allele for (allele,) in
                            execute(initialConditionsDB, "SELECT locus_allele FROM blocusalleles ORDER BY bstrain_id")]
        end
        growthrates = map(x -> fitness(x, centerallele, maxallele,
                scale, maxfitness, death, evofunction), locusalleles)
        total_abundance = sum(abundance)
        spacers = [
            Vector([spacer_id for (spacer_id,) in
                    execute(
                initialConditionsDB,
                "SELECT spacer_id FROM
                bspacers WHERE bstrain_id = $(id)
                ORDER BY spacer_id"
            )])
            for id in ids
        ]
        spacers[spacers.==[[0]]] = repeat([[]], length(spacers[spacers.==[[0]]]))
        Strains(
            maximum(ids) + 1,
            ids,
            abundance,
            total_abundance,
            spacers,
            locusalleles,
            growthrates
        )
    end
end



function make_vstrains(n_strains, n_particles_per_strain, n_pspacers_init,
    initialConditionsDB::Union{DB,Nothing})
    if initialConditionsDB == Nothing
        next_id = UInt64(n_strains + 1)
        # println("next id is $(next_id)")
        # println("type is $(typeof(next_id))")
        ids = Vector(1:n_strains)
        abundance = repeat([n_particles_per_strain], n_strains)
        total_abundance = n_strains * n_particles_per_strain
        pspacers = [
            Vector(1:n_pspacers_init) .+ repeat([n_pspacers_init * (i - 1)], n_pspacers_init)
            for i = 1:n_strains
        ]

        Strains(
            next_id,
            ids,
            abundance,
            total_abundance,
            pspacers,
            repeat([], n_strains),
            repeat([], n_strains)
        )
    else
        ids = [strain_id for (strain_id,) in
               execute(initialConditionsDB, "SELECT vstrain_id FROM vabundance ORDER BY vstrain_id")]
        abundance = [abund for (abund,) in
                     execute(initialConditionsDB, "SELECT abundance FROM vabundance ORDER BY vstrain_id")]
        total_abundance = sum(abundance)
        pspacers = [
            Vector([spacer_id for (spacer_id,) in
                    execute(
                initialConditionsDB,
                "SELECT spacer_id FROM
                vpspacers WHERE vstrain_id = $(id)
                ORDER BY spacer_id"
            )])
            for id in ids
        ]
        Strains(
            maximum(ids) + 1,
            ids,
            abundance,
            total_abundance,
            pspacers,
            repeat([], length(ids)),
            repeat([], length(ids))
        )
    end
end


#########################################################
##################THIS SAVES NEW STRAINS#################

mutable struct State
    bstrains::Strains
    vstrains::Strains
    next_pspacer_id::UInt64

    function State(
        rng, n_bstrains, n_hosts_per_bstrain, random_initial_alleles,
        initial_locus_allele, evofunction,
        maxallele, centerallele, scale, maxfitness, death,
        n_vstrains, n_particles_per_vstrain, n_pspacers_init, initialConditionsDB::Union{DB,Nothing}
    )
        if initialConditionsDB == Nothing
            next_pspacer_id = 1 + n_vstrains * n_pspacers_init
        else
            pspacers = [spacer_id for (spacer_id,) in
                        execute(initialConditionsDB,
                "SELECT DISTINCT spacer_id FROM vpspacers")]
            next_pspacer_id = UInt64(1 + maximum(pspacers))
        end
        new(
            make_bstrains(rng, n_bstrains, n_hosts_per_bstrain, random_initial_alleles, initial_locus_allele,
                evofunction, maxallele, centerallele, scale, maxfitness, death, initialConditionsDB),
            make_vstrains(n_vstrains, n_particles_per_vstrain, n_pspacers_init, initialConditionsDB),
            next_pspacer_id
        )
    end
end

function State(p::Params, initialConditionsDB::Union{DB,Nothing})
    State(
        MersenneTwister(p.rng_seed), UInt64(p.n_bstrains), UInt64(p.n_hosts_per_bstrain), 
        Bool(p.random_initial_alleles), Float64(p.initial_locus_allele),
        p.evofunction, p.max_allele, p.center_allele, p.evofunctionScale, p.max_fitness, p.microbe_death_rate,
        UInt64(p.n_vstrains), UInt64(p.n_particles_per_vstrain), UInt64(p.n_protospacers), initialConditionsDB
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

    function Simulation(p::Params, initialConditionsDB::Union{DB,Nothing})
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
        state = State(p, initialConditionsDB)
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
