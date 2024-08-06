#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]
t0 = map(x->parse(Float64,x),ARGS[2:2:end])
dt = map(x -> parse(Float64, x), ARGS[3:2:end])
##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbOutputPath = joinpath("trait-evolution_output.sqlite") # cluster
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/11_MOI3/sweep_db_gathered.sqlite") # local
# dbOutputPath = joinpath("/Volumes/Yadgah/trait-evolution.sqlite") # local

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
# execute(dbOutput, "CREATE TABLE b_fitness_richness_shift (t0 REAL, initial_richness REAL, time_elapsed REAL, final_richness REAL, richness_shift REAL)")
# execute(dbOutput, "CREATE TABLE b_fitness_shannon_shift (t0 REAL, initial_shannon REAL, time_elapsed REAL, final_shannon REAL, shannon_shift REAL)")
# execute(dbOutput, "CREATE TABLE b_fitness_expectation_shift (t0 REAL, initial_exp REAL, time_elapsed REAL, final_exp REAL, exp_shift REAL)")
# execute(dbOutput, "CREATE TABLE b_fitness_variance_shift (t0 REAL, initial_var REAL, time_elapsed REAL, final_var REAL, var_shift REAL)")


function traitEvolution()
    growthrates = DataFrame(execute(dbSim, "SELECT bstrain_id, growth_rate FROM bgrowthrates WHERE run_id = $(run_id)"))
    bfreq = innerjoin(growthrates,DataFrame(execute(dbSim, "SELECT t, bstrain_id, abundance 
            FROM babundance WHERE run_id = $(run_id)")),on=:bstrain_id)
    bfreq = innerjoin(bfreq,DataFrame(execute(dbSim, "SELECT t, microbial_abundance 
            FROM summary WHERE run_id = $(run_id)")),on=:t)
    bfreq[!, :bfreq] = bfreq[!, :abundance] ./ bfreq[!, :microbial_abundance]
    bfreq[!, :expectation] = bfreq[!, :bfreq] .* bfreq[!, :growth_rate]
    bfreq = groupby(bfreq, [:t,:growth_rate,:expectation])
    bfreq = combine(bfreq, [:bfreq .=> sum,]; renamecols=false)
    richness = groupby(bfreq, :t)
    richness = combine(richness, :growth_rate .=> length; renamecols=false)
    rename!(richness, :growth_rate => :richness)
    shannon = copy(bfreq)
    shannon[!, :shannon] = map(x->exp(x),
                -1*shannon[!, :bfreq].* map(x -> log(x), shannon[!, :bfreq]))
    shannon = groupby(shannon, :t)
    shannon = combine(shannon, [:shannon .=> prod,]; renamecols=false)
    simpson = copy(bfreq)
    simpson[!, :simpson] = simpson[!, :bfreq].^2
    simpson = groupby(simpson, :t)
    simpson = combine(simpson, [:simpson .=> sum,]; renamecols=false)
    expectation = groupby(bfreq, :t)
    expectation = combine(expectation, [:expectation .=> sum,]; renamecols=false)
    var = innerjoin(bfreq[!,[:t,:growth_rate,:bfreq]],expectation,on=:t)
    var[!, :variance] = var[!, :bfreq].*(var[!, :expectation] .- var[!, :growth_rate]) .^ 2
    var = groupby(var, :t)
    var = combine(var, [:variance .=> sum,]; renamecols=false)
    richness |> SQLite.load!(dbOutput,
        "b_intrinsic_fitness_richness", ifnotexists=true)
    shannon |> SQLite.load!(dbOutput,
        "b_intrinsic_fitness_shannon", ifnotexists=true)
    simpson |> SQLite.load!(dbOutput,
        "b_intrinsic_fitness_simpson", ifnotexists=true)
    expectation |> SQLite.load!(dbOutput,
        "b_intrinsic_fitness_expectation", ifnotexists=true)
    var |> SQLite.load!(dbOutput,
        "b_intrinsic_fitness_variance", ifnotexists=true)
    # for (k,l) in zip(t0,dt)
    #     if k+l > maximum(bfreq[!,:t]) 
    #         continue
    #     end
    #     times = unique(bfreq[!,:t])
    #     if !in(k,times)
    #         i = maximum(times[times.<k])
    #     else
    #         i = k
    #     end
    #     if !in(k+l, times)
    #         j = maximum(times[times.<k+l])
    #     else
    #         j = k+l
    #     end
    #     rich0 = richness[richness.t.==i, :].richness[1]
    #     richf = richness[richness.t.==j, :].richness[1]
    #     dRich = richf - rich0
    #     shan0 = shannon[shannon.t.==i, :].shannon[1]
    #     shanf = shannon[shannon.t.==j, :].shannon[1]
    #     dShan = shanf - shan0
    #     exp0 = expectation[expectation.t.==i, :].expectation[1]
    #     expf = expectation[expectation.t.==j, :].expectation[1]
    #     dExp = expf - exp0
    #     var0 = var[var.t.==i, :].variance[1]
    #     varf = var[var.t.==j, :].variance[1]
    #     dVar = varf - var0
    #     execute(dbOutput, "INSERT INTO b_fitness_richness_shift VALUES (?,?,?,?,?)",
    #         (k, rich0, k+l, richf, dRich))
    #     execute(dbOutput, "INSERT INTO b_fitness_shannon_shift VALUES (?,?,?,?,?)",
    #         (k, shan0, k+l, shanf, dShan))
    #     execute(dbOutput, "INSERT INTO b_fitness_expectation_shift VALUES (?,?,?,?,?)",
    #         (k, exp0, k+l, expf, dExp))
    #     execute(dbOutput, "INSERT INTO b_fitness_variance_shift VALUES (?,?,?,?,?)",
    #         (k, var0, k+l, varf, dVar))
    # end
end

traitEvolution()
println("Complete!")
