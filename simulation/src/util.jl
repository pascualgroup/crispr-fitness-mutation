
### RNG UTILITY FUNCTIONS ##

"""
Returns an index from 1:length(w) with probability proportional to w.

s must be precomputed to be sum(w).
"""
function sample_linear_integer_weights(rng::MersenneTwister, w::Vector{UInt64}, s::UInt64)
    i = rand(rng, 1:s)
    cs = 0
    for j = 1:(length(w) - 1)
        cs += w[j]
        if i <= cs
            return j
        end
    end
    length(w)
end

"""
Removes an item in the middle of an array that does not need to be kept ordered in constant time.

The item is replaced with the item at the end of the array, and then the item at the end of the
array is removed.
"""
function swap_with_end_and_remove!(a, index)
    if index != lastindex(a)
        setindex!(a, a[lastindex(a)], index)
    end
    pop!(a)
    nothing
end

function remove_strain!(strains, index)
    # This is only used when a strain has gone extinct
    @assert strains.abundance[index] == 0

    @debug "Removing strain" id=strains.ids[index] index=index

    swap_with_end_and_remove!(strains.ids, index)
    swap_with_end_and_remove!(strains.abundance, index)
    swap_with_end_and_remove!(strains.spacers, index)
    if length(strains.growthrates) > 0
        swap_with_end_and_remove!(strains.growthrates, index)
    end
    if length(strains.growthalleles) > 0
        swap_with_end_and_remove!(strains.growthalleles, index)
    end
end
