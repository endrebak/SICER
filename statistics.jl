

using Distributions: Poisson, pdf, ccdf
using StatsBase: tiedrank

include("config.jl")


function compute_window_score(i, poisson)

    if i < poisson.λ
        return 0
    end

    p_value = pdf(poisson, i)

    if p_value > 0
        window_score = -log(p_value)
    else
        # log of zero not defined
        window_score = 1000
    end

    return window_score

end


function update_island_expectations!(island_expectations, scaled_score, bin_size, poisson, island_enriched_threshold, gap_contribution)

  current_max_scaled_score = length(island_expectations) - 1
  # println("current_max_scaled_score ", current_max_scaled_score)
  # println("scaled_score ", scaled_score)
  if scaled_score > current_max_scaled_score
      #index is the scaled_score
      for index in (current_max_scaled_score + 1):(scaled_score)
          temp=0.0
          #i is the number of tags in the added windo w
          i = island_enriched_threshold
          # print("bin_size ", bin_size)
          current_island = floor(Integer, round(index - compute_window_score(i, poisson)/BIN_SIZE)) + 1
          while (current_island >= 1)
              # println("current island ", current_island)
              # println("temp ", temp)
              # println("i ", i)
              island_expectation = island_expectations[current_island]
              # if island_expectation != 0
              #     println("island_expectation: ", island_expectation, ", index ", index)
              # end

              temp += pdf(poisson, i) * island_expectation
              i += 1
              if i == 500
                  break
              end
              current_island = floor(Integer, round(index - compute_window_score(i, poisson)/BIN_SIZE)) + 1
          end
          temp *= gap_contribution
          append!(island_expectations, temp)
      end
  end

  return island_expectations

end


function generate_cumulative_distribution(island_expectations)

    l = length(island_expectations)
    cumulative = zeros(l)
    partial_sum = 0

    for i in 2:(l - 1)

        compliment = l - i
        partial_sum += island_expectations[compliment]
        cumulative[compliment] = partial_sum

    end

    partial_sum += island_expectations[l]
    cumulative[1] = partial_sum

    return cumulative

end

function compute_score_threshold(average_window_readcount,
                                 island_enriched_threshold,
                                 gap_contribution, boundary_contribution,
                                 genome_length_in_bins, bin_size)

    poisson = Poisson(average_window_readcount)

    required_p_value = pdf(poisson, island_enriched_threshold)
    # println("required_p_value ", required_p_value)

    prob = boundary_contribution * required_p_value
    # println("prob", prob)

    score = -log(required_p_value)
    # println("score ", score, " equal_to_prob")

    # println("round(score / BIN_SIZE) ", round(score / BIN_SIZE))
    current_scaled_score = trunc(Integer, round(score / BIN_SIZE))
    # println("current_scaled_score", current_scaled_score)

    island_expectations = zeros(current_scaled_score + 1)
    # println("create array")

    # println("genome_length_in_bins", genome_length_in_bins)
    island_expectations[1] = boundary_contribution * genome_length_in_bins / gap_contribution

    # println("set first")

    island_expectations[current_scaled_score + 1] = prob * genome_length_in_bins

    # println(island_expectations)
    # println("length: ", length(island_expectations))

    # println("set last")

    current_max_scaled_score = current_scaled_score

    interval = trunc(Integer, 1 / BIN_SIZE)
    partial_cumu = 0.0

    while (partial_cumu > E_VALUE_THRESHOLD || partial_cumu < 1e-100)

        # # println("Outer while")
        current_scaled_score += interval
        # println("current scaled score ", current_scaled_score)
        update_island_expectations!(island_expectations, current_scaled_score, bin_size, poisson, island_enriched_threshold, gap_contribution)
        current_expectation = island_expectations[trunc(Integer, current_scaled_score)]
        # # println("current_expectation", current_expectation)

        # # println("partial_cumu ", partial_cumu)
        l = length(island_expectations)
        if l > interval
            # # println("in if")
            # # println(island_expectations)
            # open("$(length(island_expectations))", "w+") do h
            #     write(h, string(island_expectations))
            # end

            # println("length ", length(island_expectations[(l - interval) + 1:(l - 1)]))
            # println("ie sum ", sum(island_expectations[(l - interval) + 1:(l - 1)]))
            partial_cumu = sum(island_expectations[(l - interval) + 1:(l - 1)])
        else
            partial_cumu = sum(island_expectations)
        end

    end

    cumulative = generate_cumulative_distribution(island_expectations)

    score_threshold = 0
    for (i, e) in enumerate(cumulative)
        # println("(i, e) ", (i, e))
        if e < E_VALUE
            # println("e smaller than E_VALUE. e: ", e, " i: ", i, " E_VALUE_THRESHOLD ", E_VALUE)
            score_threshold = (i - 1) * BIN_SIZE
            break
        end
    end

    return score_threshold

end


function compute_enriched_threshold(poisson)

    current_threshold, survival_function = 0, 1

    println("poisson.λ ", poisson.λ)
    survival_function -= pdf(poisson, current_threshold)
    println("survival_function ", survival_function)
    println("p ", WINDOW_P_VALUE)
    while survival_function > WINDOW_P_VALUE
        current_threshold += 1
        survival_function -= pdf(poisson, current_threshold)
    end

    island_enriched_threshold = current_threshold + 1

    return island_enriched_threshold

end


function single_gap_factor(island_enriched_threshold, poisson)

    _sum = 0
    for i in 0:(island_enriched_threshold - 1)
        # println("i ", i, " poisson: ", pdf(poisson, i))
        _sum += pdf(poisson, i)
    end

    return _sum
end


function compute_gap_factor(island_enriched_threshold, gap_intervals_allowed, poisson)

    if gap_intervals_allowed == 0
        return 1
    end


    gap_factor = single_gap_factor(island_enriched_threshold, poisson)
    # println("single_gap_factor ", gap_factor)

    max_gap_score = 1
    for i in 1:gap_intervals_allowed
        # println("i ", i, " poisson: ", gap_factor ^ i)
        max_gap_score += gap_factor ^ i
    end

    # println("max_gap_score ", max_gap_score)
    return max_gap_score
end


function compute_boundary(island_enriched_threshold, gap_intervals_allowed, poisson)

    single_gap = single_gap_factor(island_enriched_threshold, poisson)
    single_boundary_score = single_gap ^ (gap_intervals_allowed + 1)
    start_and_end_score = single_boundary_score * single_boundary_score

end



function compute_background_probabilities(total_chip_count, bin_size, effective_genome_fraction, gap_intervals_allowed)

    # println("effective_genome_fraction ", effective_genome_fraction)
    # println("total_chip_count ", total_chip_count)
    average_window_readcount = total_chip_count * (bin_size / effective_genome_fraction)
    tag_density = total_chip_count / effective_genome_fraction
    # println("tag_density ", tag_density)
    # println("tag_density * windowSize", tag_density * bin_size)


    # self.tag_density = total_tags * 1.0 / genomeLength;
    # self.average = self.tag_density * windowSize;

    # println("average_window_readcount (poisson param)", average_window_readcount)

    poisson = Poisson(average_window_readcount)

    println("compute island enriched threshold")
    island_enriched_threshold = compute_enriched_threshold(poisson)
    println("island_enriched_threshold ", island_enriched_threshold)

    println("compute gap factor")
    gap_contribution = compute_gap_factor(island_enriched_threshold, gap_intervals_allowed, poisson)
    println("gap_contribution", gap_contribution)

    println("compute boundary")
    boundary_contribution = compute_boundary(island_enriched_threshold, gap_intervals_allowed, poisson)
    println("boundary_contribution", boundary_contribution)

    genome_length_in_bins = effective_genome_fraction / bin_size

    println("compute score threshold")
    println("average_window_readcount, island_enriched_threshold, gap_contribution, boundary_contribution, genome_length_in_bins, bin_size")
    println(join([average_window_readcount, island_enriched_threshold, gap_contribution, boundary_contribution, genome_length_in_bins, bin_size], " "))
    score_threshold = compute_score_threshold(average_window_readcount, island_enriched_threshold,
                                              gap_contribution, boundary_contribution, genome_length_in_bins, bin_size)

    # println("island_enriched_threshold ", island_enriched_threshold)
    # println("gap_contribution ", gap_contribution)
    # println("boundary_contribution ", boundary_contribution)
    # println("average_window_readcount ", average_window_readcount)
    # println("score_threshold ", score_threshold)
    # println("window average", total_chip_count / effective_genome_fraction)

    return score_threshold, island_enriched_threshold, average_window_readcount

end


function give_bins_pvalues(df, island_enriched_threshold, average_window_readcount)

    println("df give bins p", head(df))
    println("island_enriched_threshold", island_enriched_threshold)
    df = df[df[:Count] .>= island_enriched_threshold, :]
    println("df give bins p", head(df))

    println("poisson value give bins p-values ", average_window_readcount)
    poisson = Poisson(average_window_readcount)
    df[:Score] = pdf.(poisson, df[:Count])
    println("scores: ", unique(df[:Score]))
    println("counts: ", unique(df[:Count]))
    df
end


function give_islands_fdr_score(df, total_chip_reads, total_input_reads, effective_genome_fraction)

    # println("total_chip_reads ", total_chip_reads)
    # println("total_input_reads ", total_input_reads)
    # total_chip_reads = 9924

    scaling_factor = (total_chip_reads) / total_input_reads
    zero_controls_multiplier = total_input_reads / effective_genome_fraction

    avg_0_denom = (df[:End] - df[:Start] .+ 1) .* zero_controls_multiplier
    println("avg_0_denom")
    println(typeof(avg_0_denom))
    println(avg_0_denom[1:5])
    avg_0_denom[avg_0_denom .> 0.25] .= 0.25
    println(avg_0_denom[1:5])
    avg_0_denom = avg_0_denom .* scaling_factor
    println(avg_0_denom[1:5])

    avg = df.InputCount .* scaling_factor
    println("avg_0_denom[df.InputCount .== 0][1:5]", avg_0_denom[df.InputCount .== 0][1:5])
    println("avg[df.InputCount .== 0][1:5]", avg[df.InputCount .== 0][1:5][1:5])
    avg[df.InputCount .== 0] = avg_0_denom[df.InputCount .== 0]
    println("avg_0_denom[df.InputCount .== 0][1:5]", avg_0_denom[df.InputCount .== 0][1:5])
    println("avg[df.InputCount .== 0][1:5]", avg[df.InputCount .== 0][1:5][1:5])

    poisson = Poisson.(avg)

    fold_change = log2.(df[:Count] ./ avg)

    # println("fold_change")
    # println(fold_change)
    df[:log2FoldChange] = fold_change

    p_values = DataFrame(PValue = ccdf.(poisson, df.Count))

    println("p_values")
    # println(p_values)
    println(typeof(p_values))
    # DataFrame

    println("avg")
    println(typeof(avg))
    println(length(avg))
    # println(df.Count)
    println("length(avg) ", length(avg))
    println("typeof(df.Count <= avg) ", typeof(df.Count <= avg))
    # println("length(typeof(df.Count <= avg)) ", length(typeof(df.Count <= avg)))
    println("typeof(df.Count .<= avg) ", typeof(df.Count .<= avg))
    println("df.Count .<= avg) ", df.Count .<= avg)
    # println("length(df.Count .<= avg) ", typeof(df.Count .<= avg))
    boolean_index = df.Count .<= avg
    println("len bool idx", length(boolean_index))
    println("len p_values", length(p_values))
    p_values[boolean_index, 1] .= 1

    println(head(p_values))

    println(head(df))
    df[:PValue] = p_values[1]
    println(head(df))

    fdr = (p_values[1] .* nrow(df)) ./ tiedrank(p_values[1])
    fdr[fdr .> 1, :] .= 1

    println("fdr ", fdr[1:10])
    df[:FDR] = fdr

    println("head(df) ", head(df))

    return df
end
