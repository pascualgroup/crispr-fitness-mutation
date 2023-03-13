#!/usr/bin/env julia

function currentStructure!(matchStructure::tripartite)
    time = matchStructure.time
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    for (vstrain_id,) in execute(dbTempMatch,"SELECT DISTINCT vstrain_id
            FROM vstrain_matched_pspacers WHERE t = $(time)
            ORDER BY vstrain_id")
        matchStructure.strainID = Int64(vstrain_id)
        phenotype = [type for (type,) in
        execute(dbTempMatch,"SELECT matched_pspacer_id
        FROM vstrain_matched_pspacers
        WHERE t = $(time) AND vstrain_id = $(vstrain_id)
        ORDER BY matched_pspacer_id")]
        matchStructure.phenotype = phenotype
        # println("phenotype is $(phenotype)")
        if in(phenotype,matchStructure.vphenotypes)
            if in(phenotype,keys(matchStructure.vstrainclasses))
                push!(matchStructure.vstrainclasses[phenotype],Int64(vstrain_id))
            else
                newCurrent!(matchStructure)
            end
        else
            newPhenotype!(matchStructure)
        end

    end
    matchStructure.vb = false
    for (bstrain_id,) in execute(dbTempMatch,"SELECT DISTINCT bstrain_id
            FROM bstrain_matched_spacers WHERE t = $(time)
            ORDER BY bstrain_id")
        if length([strain for (strain,) in
                execute(dbTempMatch,"SELECT vstrain_id
                    FROM bstrain_to_vstrain_matches
                    WHERE t = $(time) AND match_length = 1
                    AND bstrain_id = $(bstrain_id)")]) == 0
            continue
        end
        matchStructure.strainID = Int64(bstrain_id)
        phenotype = [type for (type,) in
        execute(dbTempMatch,"SELECT matched_spacer_id
        FROM bstrain_matched_spacers
        WHERE t = $(time) AND bstrain_id = $(bstrain_id) 
        ORDER BY matched_spacer_id")]
        matchStructure.phenotype = phenotype
        # println("phenotype is $(phenotype)")
        if in(phenotype,matchStructure.bphenotypes)
            if in(phenotype,keys(matchStructure.bstrainclasses))
                push!(matchStructure.bstrainclasses[phenotype],Int64(bstrain_id))
            else
                newCurrent!(matchStructure)
            end
        else
            newPhenotype!(matchStructure)
        end
    end
end

function identifyCurrentMatches!(matchStructure::tripartite)
    vb = matchStructure.vb
    time = matchStructure.time
    dbTempSim = matchStructure.dbSim
    dbTempMatch = matchStructure.dbMatch
    phenotype = matchStructure.phenotype
    strain_id = matchStructure.strainID
    if vb
        matchStructure.sbstrains[phenotype] = [Int64(strainID)
            for (strainID,) in execute(dbTempMatch,
                "SELECT bstrain_id FROM bstrain_to_vstrain_0matches
                WHERE t = $(time) AND vstrain_id = $(strain_id) 
                ORDER BY bstrain_id" )]
        if length(matchStructure.sbstrains[phenotype]) == 0
            matchStructure.sBiomass[phenotype] = 0
        else
            matchStructure.sBiomass[phenotype] = sum([Int64(abund)
                for (abund,) in execute(dbTempSim,"SELECT abundance
                    FROM babundance WHERE t = $(time) AND bstrain_id in
                    ($(join(matchStructure.sbstrains[phenotype],", ")))")])
        end
        allBstrains = [Int64(strain)for (strain,) in execute(dbTempSim,
                        "SELECT bstrain_id FROM babundance 
                        WHERE t = $(time)") if strain != 1]
        matchStructure.ibstrains[phenotype] =
            setdiff(allBstrains, matchStructure.sbstrains[phenotype])
        if length(matchStructure.ibstrains[phenotype]) == 0
            matchStructure.iBiomass[phenotype] = 0
        else
            matchStructure.iBiomass[phenotype] = sum([Int64(abund)
                for (abund,) in execute(dbTempSim, "SELECT abundance
                    FROM babundance WHERE t = $(time) AND bstrain_id in
                    ($(join(matchStructure.ibstrains[phenotype],", ")))")])
        end
    end
end