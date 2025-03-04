#!/usr/bin/env julia

println("(Julia compilation delay...)")

using SQLite
using DataFrames
using SQLite.DBInterface: execute

run_id = ARGS[1]
uPercents = []
lPercents = []
for i in ARGS[2:end-1]
    if i[1] == 'u'
        append!(uPercents, parse(Float64,chop(i, head=1,tail=0))/100)
    end
    if i[1] == 'l'
        append!(lPercents, parse(Float64, chop(i, head=1,tail=0))/100)
    end
end
# uPercent = [chop(ARGS[2],head=1) for i in ARGS[3:end-1]]
# uPercent = parse(Float64, ARGS[2]) / 100
# # lPercent = parse(Float64, ARGS[3]) / 100
# lPercents = [parse(Float64, i) / 100 for i in ARGS[3:end-1]]
shannonThreshold = parse(Float64, ARGS[end])

##
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))

dbOutputPath = joinpath("walls-shannon_output.sqlite") # cluster
# dbOutputPath = joinpath("/Volumes/Yadgah/avalanches_output.sqlite") # local
rm(dbOutputPath, force=true)
if isfile(dbOutputPath)
    error("$(dbOutputPath) already exists")
end

dbSimPath = joinpath(SCRIPT_PATH, "..", "..", "simulation", "sweep_db_gathered.sqlite") # cluster
# dbSimPath = joinpath("/Volumes/Yadgah/sylvain-martin-collab/21_MOI3", "sweep_db_gathered.sqlite")
# dbSimPath = joinpath("/Volumes/Yadgah/avalanches/simulations/comboID-1.sqlite")

dbSim = SQLite.DB(dbSimPath)
dbOutput = SQLite.DB(dbOutputPath)
execute(dbOutput, "CREATE TABLE shannon_diversity (t REAL, vshannon REAL, bshannon REAL)")
execute(dbOutput, "CREATE TABLE richness (t REAL, vrichness REAL, brichness REAL)")

# (cID,) = execute(dbSim, "SELECT combo_id FROM runs WHERE run_id = $(run_id)")
# cID = cID.combo_id
# (replicate,) = execute(dbSim, "SELECT replicate FROM runs WHERE run_id = $(run_id)")
# replicate = replicate.replicate

tmpPath = joinpath(ENV["SLURM_TMPDIR"], "walls-$(run_id).sqlite")
# tmpPath = joinpath("/Volumes/Yadgah/test.sqlite") # local
println("this is the temp path: $(tmpPath)")
rm(tmpPath, force=true)
dbTemp = SQLite.DB(tmpPath)
# dbTemp = SQLite.DB()
execute(dbTemp, "CREATE TABLE summary (t REAL, microbial_abundance INTEGER, viral_abundance INTEGER)")
execute(dbTemp, "CREATE TABLE babundance (t REAL, bstrain_id INTEGER, abundance INTEGER)")
execute(dbTemp, "CREATE TABLE vabundance (t REAL, vstrain_id INTEGER, abundance INTEGER)")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "ATTACH DATABASE '$(dbSimPath)' as dbSim")
execute(dbTemp, "INSERT INTO summary(t, microbial_abundance,viral_abundance) SELECT t, microbial_abundance, viral_abundance FROM dbSim.summary WHERE run_id = $(run_id);")
execute(dbTemp, "INSERT INTO babundance (t, bstrain_id, abundance) SELECT t, bstrain_id, abundance FROM dbSim.babundance WHERE run_id = $(run_id);")
execute(dbTemp, "INSERT INTO vabundance (t, vstrain_id, abundance) SELECT t, vstrain_id, abundance FROM dbSim.vabundance WHERE run_id = $(run_id);")
execute(dbTemp, "COMMIT")

execute(dbTemp, "BEGIN TRANSACTION")
execute(dbTemp, "CREATE INDEX summary_index ON summary (t,microbial_abundance,viral_abundance)")
execute(dbTemp, "CREATE INDEX bstrain_index ON babundance (t,bstrain_id)")
execute(dbTemp, "CREATE INDEX vstrain_index ON vabundance (t,vstrain_id)")
execute(dbTemp, "COMMIT")
# println("Processing shannon entropies of run $(run_id)")


function shannon(dbOutput, dbTemp)
    execute(dbOutput, "BEGIN TRANSACTION")
    for (time,) in execute(dbTemp, "SELECT DISTINCT t FROM summary")
        (totAbunds,) = execute(dbTemp, "SELECT microbial_abundance,viral_abundance FROM summary WHERE t = ?", (time,))
        btotal = totAbunds.microbial_abundance
        vtotal = totAbunds.viral_abundance

        bsubAbunds = [bsubAbunds.abundance for bsubAbunds in execute(dbTemp, "SELECT abundance FROM babundance WHERE t = ?", (time,))]
        brelAbunds = bsubAbunds ./ btotal
        brelAbunds = brelAbunds[brelAbunds.>=1e-200]

        if btotal == 0
            bshannon = 0
            brichness = 0
        else
            bshannon = exp(-1 * sum(brelAbunds .* (log.(brelAbunds))))
            brichness = length(brelAbunds)
        end

        vsubAbunds = [vsubAbunds.abundance for vsubAbunds in execute(dbTemp, "SELECT abundance FROM vabundance WHERE t = ?", (time,))]
        vrelAbunds = vsubAbunds ./ vtotal
        vrelAbunds = vrelAbunds[vrelAbunds.>=1e-200]

        if vtotal == 0
            vshannon = 0
            vrichness = 0
        else
            vshannon = exp(-1 * sum(vrelAbunds .* (log.(vrelAbunds))))
            vrichness = length(vrelAbunds)
        end

        execute(dbOutput, "INSERT INTO shannon_diversity VALUES (?,?,?)", (time, vshannon, bshannon))
        execute(dbOutput, "INSERT INTO richness VALUES (?,?,?)", (time, vrichness, brichness))
    end
    execute(dbOutput, "COMMIT")
end

# This function counts number of peaks, logs time series of peaks their respective durations
function peakwallCount(uPercent, lPercent, shannonThreshold, dbOutput, dbSim, dbTemp)
    (comboID,) = execute(dbSim, "SELECT combo_id FROM runs WHERE run_id = $(run_id)")
    comboID = comboID.combo_id
    (cCapacity,) = execute(dbSim, "SELECT microbe_carrying_capacity FROM param_combos WHERE combo_id = $(comboID)")
    cCapacity = cCapacity.microbe_carrying_capacity

    # Save as data frame
    series = DataFrame(execute(dbTemp, "SELECT t, microbial_abundance FROM summary"))

    upperThreshold = floor(cCapacity * uPercent)
    lowerThreshold = floor(cCapacity * lPercent)

    thresholdVals = DataFrame(carrying_capacity=cCapacity,
        upper_percent=uPercent, lower_percent=lPercent,
        upper_threshold=upperThreshold, lower_threshold=lowerThreshold,
        shannon_threshold=shannonThreshold)
    thresholdVals |> SQLite.load!(dbOutput, "threshold_values", ifnotexists=true)

    # Identify crossings of thresholds
    coarseSeries = series[(series.microbial_abundance.>=upperThreshold).|(series.microbial_abundance.<=lowerThreshold), :]
    upper = 2 * (series.microbial_abundance .>= upperThreshold) # 2 signifies values above upper threshold
    lower = series.microbial_abundance .<= lowerThreshold # 1 signifies values below lower threshold, 0 for values between lower and upper threshold,
    coarseSeries = upper + lower

    # Initialize each table of output database as data frame
    peakSeries = DataFrame(t=series.t, peak_presence=Array{Int64,1}(zeros(length(series.t))))
    peakDurations = DataFrame(begin_row=Int64[], end_row=Int64[], begin_t=Float64[], end_t=Float64[],
        peak_number=Int64[], wall_number=Int64[], duration=Float64[])
    wallSeries = DataFrame(t=series.t, wall_presence=Array{Int64,1}(zeros(length(series.t))), bshannon=Array{Float64,1}(zeros(length(series.t))), vshannon=Array{Float64,1}(zeros(length(series.t))))

    peaks = 0
    walls = 0
    check = 0
    j1, j2, j3, j4 = 0, 0, 0, 0

    for i in collect(1:length(coarseSeries))
        println("Scanning through time series point $(i)")
        if i == length(coarseSeries)
            return peaks, walls, peakSeries, peakDurations, wallSeries
        end
        if coarseSeries[i] == 1 && coarseSeries[i+1] == 0
            check = 1
            j1 = i
        end
        if coarseSeries[i] == 0 && coarseSeries[i+1] == 2
            if check == 1
                check = 2
                j2 = i
            end
        end
        if coarseSeries[i] == 2 && coarseSeries[i+1] == 0
            if check == 2
                check = 3
                j3 = i
            end
        end
        if coarseSeries[i] == 0 && coarseSeries[i+1] == 1
            if check == 3
                j4 = i
                peaks = peaks + 1

                peakSeries.peak_presence[j1:j4] = Array{Int64,1}(ones(length(j1:j4)))
                peak_duration = peakSeries.t[j4] - peakSeries.t[j1]

                timeStmt = join(map(x -> string(" $(x)"), peakSeries.t[(j1):j4]), ',')
                # timeStmt = append!([string("t = $(peakSeries.t[j1])")], timeStmt)
                bshannon = [shannon.bshannon for shannon in execute(dbOutput, "SELECT bshannon FROM shannon_diversity WHERE t in ($(timeStmt...)) ORDER BY t")]
                vshannon = [shannon.vshannon for shannon in execute(dbOutput, "SELECT vshannon FROM shannon_diversity WHERE t in ($(timeStmt...)) ORDER BY t")]

                if maximum(bshannon) >= shannonThreshold
                    walls += 1
                    wallSeries.wall_presence[j1:j4] = Array{Int64,1}(ones(length(j1:j4)))
                    wallSeries.bshannon[j1:j4] = bshannon[:, :]
                    wallSeries.vshannon[j1:j4] = vshannon[:, :]
                    push!(peakDurations, [j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks walls peak_duration])
                else
                    push!(peakDurations, [j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks 0 peak_duration])
                end

                check = 0
            end
        end
    end
    return peaks, walls, peakSeries, peakDurations, wallSeries
end

function walls(uPercents, lPercents, shannonThreshold, dbOutput, dbSim, dbTemp)
    shannon(dbOutput, dbTemp)
    for uPercent in uPercents 
        for lPercent in lPercents 
            peaks, walls, peakSeries, peakDurations, wallSeries = peakwallCount(uPercent, lPercent, shannonThreshold, dbOutput, dbSim, dbTemp);
            count = DataFrame(num_peaks=peaks, num_walls=walls)
            #
            count[!, :upper_percent] .= uPercent
            peakSeries[!, :upper_percent] .= uPercent
            peakDurations[!, :upper_percent] .= uPercent
            wallSeries[!, :upper_percent] .= uPercent
            #
            count[!, :lower_percent] .= lPercent
            peakSeries[!, :lower_percent] .= lPercent
            peakDurations[!, :lower_percent] .= lPercent
            wallSeries[!, :lower_percent] .= lPercent

            # peakSeries |> SQLite.load!(dbOutput, "microbial_peak_series", ifnotexists=true)
            peakDurations |> SQLite.load!(dbOutput, "microbial_peakwall_durations", ifnotexists=true)
            count |> SQLite.load!(dbOutput, "microbial_peakwall_count", ifnotexists=true)
            # wallSeries |> SQLite.load!(dbOutput, "microbial_wall_series", ifnotexists=true)
        end
    end
end

walls(uPercents, lPercents, shannonThreshold, dbOutput, dbSim, dbTemp)

rm(tmpPath, force=true)
println("Complete!")
