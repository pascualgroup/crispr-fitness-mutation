function setMicrobeCommunity(fxnNum, auxiliary, rngDiv, n_strains, n_hosts_per_strain,
    vpspacers,
    initial_locus_allele, evofunction,
    maxallele, locusresolution, centerallele, scale, maxfitness, death)
    if fxnNum == 0
        initialalleles = repeat([initial_locus_allele], n_strains)
        initialgrowthrates = map(x -> fitness(x, centerallele, maxallele,
                scale, maxfitness, death, evofunction), initialalleles)
        ids = Vector(1:n_strains)
        abundance = repeat([n_hosts_per_strain], n_strains)
        total_abundance = n_strains * n_hosts_per_strain
        return Strains(
            n_strains + 1, # without initial conditions we just have one strain, and the next_strain_id is just 2
            ids,
            abundance,
            total_abundance,
            repeat([[]], n_strains),
            initialalleles,
            initialgrowthrates
        )
    end
    if fxnNum == 1
        minSelectedAllele = 2 * centerallele - maxallele # this should be greater than center_allele - max_allele
        maxSelectedAllele = maxallele# this should be less than max_allele
        alleles = abs.(collect(minSelectedAllele:locusresolution:maxSelectedAllele))
        # locusalleles = [centerallele, sample(rngDiv, alleles, n_strains - 1; replace=false)...] # this ensures that source of fuel has high growth rate
        locusalleles = collect(centerallele:locusresolution:maxSelectedAllele)
        rem = length(locusalleles) % n_strains
        maxIdx = length(locusalleles) - rem
        idx = collect(Int64(1):Int64(maxIdx / n_strains):Int64(maxIdx))
        locusalleles = locusalleles[idx]
        growthrates = Vector(map(x -> fitness(x, centerallele, maxallele,
                scale, maxfitness, death, evofunction), locusalleles))
        total_abundance = n_strains * n_hosts_per_strain
        freqs = repeat([1 / n_strains], n_strains)
        ids = collect(1:n_strains)
        #
        singlematches = sample(rngDiv, vpspacers..., Int64(n_strains - 1); replace=false)
        spacers = [Vector([]), map(x -> Vector([x]), singlematches)...]
        abundance = Vector(map(x -> floor(x * total_abundance), freqs))
        #
        return Strains(
            maximum(ids) + 1,
            ids,
            abundance,
            total_abundance,
            spacers,
            locusalleles,
            growthrates
        )
    end
    if fxnNum == 2
        minSelectedAllele = 2 * centerallele - maxallele # this should be greater than center_allele - max_allele
        maxSelectedAllele = maxallele# this should be less than max_allele
        alleles = abs.(collect(minSelectedAllele:locusresolution:maxSelectedAllele))
        # locusalleles = [centerallele, sample(rngDiv, alleles, n_strains - 1; replace=false)...] # this ensures that source of fuel has high growth rate
        locusalleles = collect(centerallele:locusresolution:maxSelectedAllele)
        rem = length(locusalleles) % n_strains
        maxIdx = length(locusalleles) - rem
        idx = collect(Int64(1):Int64(maxIdx / n_strains):Int64(maxIdx))
        locusalleles = locusalleles[idx]
        println("locus alleles: $(locusalleles)")
        growthrates = Vector(map(x -> fitness(x, centerallele, maxallele,
                scale, maxfitness, death, evofunction), locusalleles))
        #
        total_abundance = n_strains * n_hosts_per_strain
        freqs = repeat([1 / n_strains], n_strains)
        #
        ids = collect(1:n_strains)
        spacers = [[rand(rngDiv, vpspacers[i]) for i in 1:length(vpspacers) if i .!= j] for j in 1:n_strains]
        abundance = Vector(map(x -> floor(x * total_abundance), freqs))
        #
        println("total abundance of microbes: $(total_abundance)")
        return Strains(
            maximum(ids) + 1,
            ids,
            abundance,
            total_abundance,
            spacers,
            locusalleles,
            growthrates
        )
    end
    if fxnNum == 3
        commSampleDB = DB(auxiliary[1])
        runID = auxiliary[2]
        sampleT = auxiliary[3]
        ids = [strain_id for (strain_id,) in
                execute(commSampleDB, "SELECT bstrain_id FROM babundance 
                        ORDER BY bstrain_id 
                        WHERE run_id = $(runID) AND t = $(sampleT)")]
        abundance = [abund for (abund,) in
                        execute(commSampleDB, "SELECT abundance FROM babundance 
                        ORDER BY bstrain_id 
                        WHERE run_id = $(runID) AND t = $(sampleT)")]
        locusalleles = [allele for (allele,) in
                        execute(commSampleDB, "SELECT locus_allele FROM blocusalleles 
                        ORDER BY bstrain_id
                        WHERE run_id = $(runID) 
                        AND bstrain_id in ($(join(ids,", ")))"
        )]
        growthrates = map(x -> fitness(x, centerallele, maxallele,
                scale, maxfitness, death, evofunction), locusalleles)
        total_abundance = sum(abundance)
        spacers = [
            Vector{UInt64}([spacer_id for (spacer_id,) in
                    execute(
                commSampleDB,
                "SELECT spacer_id FROM
                bspacers 
                WHERE run_id = $(runID) 
                AND bstrain_id = $(id)
                ORDER BY spacer_id"
            )])
            for id in ids
        ]
        return Strains(
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


function setVirusCommunity(fxnNum, auxiliary, n_strains, n_particles_per_strain, n_pspacers_init)
    if fxnNum == 0
        next_id = UInt64(n_strains + 1)
        next_pspacer_id = 1 + n_strains * n_pspacers_init
        ids = Vector(1:n_strains)
        abundance = repeat([n_particles_per_strain], n_strains)
        total_abundance = n_strains * n_particles_per_strain
        pspacers = [
            Vector(1:n_pspacers_init) .+ repeat([n_pspacers_init * (i - 1)], n_pspacers_init)
            for i = 1:n_strains
        ]

        return next_pspacer_id, Strains(
            next_id,
            ids,
            abundance,
            total_abundance,
            pspacers,
            repeat([], n_strains),
            repeat([], n_strains)
        )
    end
    if fxnNum == 1
        abundance = [n_particles_per_strain]
        total_abundance = n_particles_per_strain
        ids = [1]
        pspacers = [Vector(collect(1:n_pspacers_init))]
        next_pspacer_id = UInt64(1 + maximum(pspacers[1]))
        return next_pspacer_id, Strains(
            maximum(ids) + 1,
            ids,
            abundance,
            total_abundance,
            pspacers,
            repeat([], length(ids)),
            repeat([], length(ids))
        )
    end
    if fxnNum == 2
        total_abundance = Int64(floor(n_particles_per_strain/n_strains)) * n_strains
        abundance = repeat([Int64(floor(n_particles_per_strain/n_strains))], n_strains)
        ids = collect(1:n_strains)
        pspacers = [
            Vector(1:n_pspacers_init) .+ repeat([n_pspacers_init * (i - 1)], n_pspacers_init)
            for i = 1:n_strains
        ]
        next_pspacer_id = UInt64(1 + maximum(reduce(vcat,pspacers)))
        return next_pspacer_id, Strains(
            maximum(ids) + 1,
            ids,
            abundance,
            total_abundance,
            pspacers,
            repeat([], length(ids)),
            repeat([], length(ids))
        )
    end
    if fxnNum == 3
        commSampleDB = DB(auxiliary[1])
        runID = auxiliary[2]
        sampleT = auxiliary[3]
        ids = [strain_id for (strain_id,) in
                execute(commSampleDB, "SELECT vstrain_id FROM vabundance 
                        ORDER BY vstrain_id 
                        WHERE run_id = $(runID) AND t = $(sampleT)")]
        abundance = [abund for (abund,) in
                        execute(commSampleDB, "SELECT abundance FROM vabundance 
                        ORDER BY vstrain_id 
                        WHERE run_id = $(runID) AND t = $(sampleT)")]
        total_abundance = sum(abundance)
        pspacers = [
            Vector{UInt64}([spacer_id for (spacer_id,) in
                    execute(
                commSampleDB,
                "SELECT spacer_id FROM
                vpspacers 
                WHERE run_id = $(runID) 
                AND vstrain_id = $(id)
                ORDER BY spacer_id"
            )])
            for id in ids
        ]
        next_pspacer_id = UInt64(1 + maximum(reduce(vcat, pspacers)))
        return next_pspacer_id, Strains(
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