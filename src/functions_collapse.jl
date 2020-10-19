using StatsBase
using Roots

include("functions_align.jl")

# reduced_is_same restricts comparisons to sites that are different
# from the reference in either sequence being compared.
# NB this could be changed to the consensus sequence for the alignment to avoid
# reliance on an another file
function reduced_is_same(nuc_bit_array, ref_array)
    same_array_size = size(nuc_bit_array, 2)

    same_array = trues(same_array_size, same_array_size)

    println("number of threads for sequence comparisons: ", Threads.nthreads())

    difference_list = get_list_of_differences(nuc_bit_array, ref_array)

    for a in 1:size(nuc_bit_array,2) - 1

        for b in (a + 1):size(nuc_bit_array,2)

            start_again = true

            for r in difference_list[a]
                @inbounds x = nuc_bit_array[r,a]
                @inbounds y = nuc_bit_array[r,b]
                different = (x & y) < 16

                # store the different True/False result in an array:
                if different
                    same_array[a,b] = !different
                    same_array[b,a] = !different
                    start_again = false
                    break
                end
            end

            # need to check the second seq's sites too, if we get this far:
            if start_again
                for r in difference_list[b]
                    if r in difference_list[a]
                        continue
                    end

                    @inbounds x = nuc_bit_array[r,a]
                    @inbounds y = nuc_bit_array[r,b]
                    different = (x & y) < 16

                    # store the different True/False result in an array:
                    if different
                        same_array[a,b] = !different
                        same_array[b,a] = !different
                        break
                    end
                end
            end
        end
    end

    return same_array
end

# # compare sites (rows) until a difference is found, for all pairs of sequences (columns)
# # - break as soon as you know they're different - don't need any more information from this pair
# function is_same(nuc_bit_array)
#     same_array_size = size(nuc_bit_array, 2)
#
#     same_array = trues(same_array_size, same_array_size)
#
#     println("number of threads for sequence comparisons: ", Threads.nthreads())
#
#     for a in 1:size(nuc_bit_array,2) - 1
#         for b in (a + 1):size(nuc_bit_array,2)
#             for r in Iterators.reverse(1:size(nuc_bit_array,1))
#                 @inbounds x = nuc_bit_array[r,a]
#                 @inbounds y = nuc_bit_array[r,b]
#                 different = (x & y) < 16
#
#                 # store the different True/False result in an array:
#                 if different
#                     same_array[a,b] = !different
#                     same_array[b,a] = !different
#                     break
#                 end
#             end
#         end
#     end
#
#     return same_array
# end


# Functions for parallelising
# ---------------------------

# This to provide the chunks to chop up the n*n matrix into approximately equal
# amounts of work per thread when we parallelise it.

# Because we are effectively only working on one triangle of a square matrix,
# the first rows have many more columns to fill, so that giving threads
# equal numbers of rows to populate means very unequal amounts of work, per thread.

# Instead, we chop up the triangle into equal-area trapezoids.


# For a given sample size (n) and quantile (q), generate an expression that
# we can find the root for (solve for x) - this will provide the row number
# that splits a right-angled triangle of area n^2 / 2 into quantiles q & 1-q by
# area:
function make_equation(n, q)
    f = x-> x * n  - ((x^2 + x) / 2) - ((n^2 - n) / 2) * q
    return f
end

# get many splits:
function get_chunk_splits(n, q_Array)
    chunk_splits = Array{Int, 1}(undef, length(q_Array))
    for (i, q) in enumerate(q_Array)
        f = make_equation(n, q)
        split = find_zero(f, (0, n), Bisection())
        chunk_splits[i] = floor(split)
    end

    return chunk_splits
end

# turn the splits into chunks that are ranges to iterate over:
function get_ranges(n, N_threads)
    q_Array = Array{Float64, 1}(undef, N_threads - 1)
    for i in Iterators.reverse(1:N_threads - 1)
        q_Array[i] = i/N_threads
    end

    chunk_splits = get_chunk_splits(n, q_Array)

    ranges = Array{Tuple{Int, Int}, 1}(undef, N_threads)
    previous_finish = 0
    for (i,split) in enumerate(chunk_splits)
        start = previous_finish + 1
        finish = split
        ranges[i] = (start, finish)
        previous_finish = finish
    end
    ranges[lastindex(ranges)] = (previous_finish + 1, n - 1)

    return ranges
end

function parallel_reduced_is_same(nuc_bit_array, ref_array)

    println("number of threads for sequence comparisons: ", Threads.nthreads())

    same_array_size = size(nuc_bit_array, 2)

    same_array = Array{Bool,2}(undef, same_array_size, same_array_size)

    for i in 1:size(same_array, 2)
        same_array[i,i] = true
    end

    # make as many chunks as we have threads
    n_Chunks = Threads.nthreads()

    # ...unless we have too many threads for this many sequences, in which case, down-size.
    # A chunk size less than or equal to half the total number of sequences makes valid ranges
    if n_Chunks > size(nuc_bit_array, 2) / 2
        n_Chunks::Int = floor(size(nuc_bit_array, 2) / 2)
    end

    ranges = get_ranges(size(nuc_bit_array, 2), n_Chunks)

    difference_list = get_list_of_differences(nuc_bit_array, ref_array)

    Threads.@threads for range in ranges
        for a in range[1]:range[2]
            for b in (a + 1):size(nuc_bit_array,2)

                start_again = true
                all_same = true

                for r in difference_list[a]
                    @inbounds x = nuc_bit_array[r,a]
                    @inbounds y = nuc_bit_array[r,b]
                    different = (x & y) < 16

                    # store the different True/False result in an array:
                    if different
                        all_same = false
                        start_again = false
                        break
                    end
                end

                # need to check the second seq's sites too, if we get this far:
                if start_again
                    for r in difference_list[b]
                        if r in difference_list[a]
                            continue
                        end

                        @inbounds x = nuc_bit_array[r,a]
                        @inbounds y = nuc_bit_array[r,b]
                        different = (x & y) < 16

                        # store the different True/False result in an array:
                        if different
                            all_same = false
                            break
                        end
                    end
                end

                same_array[a,b] = all_same
                same_array[b,a] = all_same
            end
        end
    end

    return same_array
end

# # A parallel version of is_same()
# function parallel_is_same(nuc_bit_array)
#     same_array_size = size(nuc_bit_array, 2)
#
#     # initialises as falses so have to change this:
#     # same_array = SharedArray{Bool,2}((same_array_size,same_array_size))
#     same_array = Array{Bool,2}(undef, same_array_size, same_array_size)
#
#     for i in 1:size(same_array, 2)
#         same_array[i,i] = true
#     end
#
#     # make as many chunks as we have threads
#     n_Chunks = Threads.nthreads()
#
#     # ...unless we have too many threads for this many sequences, in which case, down-size.
#     # A chunk size less than or equal to half the total number of sequences makes valid ranges
#     if n_Chunks > size(nuc_bit_array, 2) / 2
#         n_Chunks::Int = floor(size(nuc_bit_array, 2) / 2)
#     end
#
#     ranges = get_ranges(size(nuc_bit_array, 2), n_Chunks)
#
#     println("number of threads for sequence comparisons: ", Threads.nthreads())
#     # println("number of chunks: ", length(ranges))
#
#     Threads.@threads for range in ranges
#         for a in range[1]:range[2]
#             for b in (a + 1):size(nuc_bit_array,2)
#                 # comb_test+=1
#                 all_same = true
#                 for r in Iterators.reverse(1:size(nuc_bit_array,1))
#                     @inbounds x = nuc_bit_array[r,a]
#                     @inbounds y = nuc_bit_array[r,b]
#                     different = (x & y) < 16
#
#                     # store the different True/False result in an array:
#                     if different
#                         all_same = false
#                         break
#                     end
#                 end
#                 same_array[a,b] = all_same
#                 same_array[b,a] = all_same
#             end
#         end
#     end
#
#     return same_array
# end

# end of functions for parallelising
# ---------------------------

function get_one_set_from_view(bool_view)
    # one subset of rows of the view
    A = Array{Int64,1}()
    original_row_numbers = bool_view.indices[1]

    for (index, value) in enumerate(bool_view)
        if value
            push!(A, original_row_numbers[index])
        end
    end
    return Set(A)
end

function get_sets_from_view(bool_view)
    S = Set{Set}()
    for i in 1:size(bool_view, 2)
        push!(S, get_one_set_from_view(view(bool_view, :, i)))
    end
    return collect(S)
end

function check_the_view_is_good(a_bool_array)
    # because these logical arrays have symmetry only need
    # to check one triangle, but haven't here
    for i in a_bool_array
        if !i
            return false
        end
    end
    return true
end

# function is_subset(s1, s2)
#     test = true
#     for x in s1
#         if ! (x in s2)
#             test = false
#             break
#         end
#     end
#     return test
# end

# # test all sets for being supersets, and throw them out if they are
# # NB - this is 1/10 of the run time and allocates a lot of memory in total
# function get_subsets(array_of_sets)
#
#     A = Array{BitSet,1}()
#     sizehint!(A, length(array_of_sets))
#
#     for i in 1:size(array_of_sets, 1)
#         valid = true
#         for j in 1:size(array_of_sets, 1)
#             if i == j
#                 continue
#             end
#             # if array_of_sets[i] is a superset of comparison set, then break, else push to A
#             # (note the reverse logic here, because there is no superset function)
#             if issubset(array_of_sets[j], array_of_sets[i]) && array_of_sets[i] != array_of_sets[j]
#                 valid = false
#                 break
#             end
#         end
#
#         if valid
#             push!(A, array_of_sets[i])
#         end
#     end
#
#     return A
# end

# a parallel version of get_subsets()
function parallel_get_subsets(array_of_sets)

    A = Array{BitSet,1}()
    sizehint!(A, length(array_of_sets))

    tests = Array{Bool, 1}(undef, length(array_of_sets))

    for i in eachindex(tests)
        tests[i] = true
    end

    Threads.@threads for i in 1:size(array_of_sets, 1)
        for j in 1:size(array_of_sets, 1)
            if i == j
                continue
            end
            # if array_of_sets[i] is a superset of comparison set, then break, else push to A
            # (note the reverse logic here, because there is no superset function)
            if issubset(array_of_sets[j], array_of_sets[i]) && array_of_sets[i] != array_of_sets[j]
                tests[i] = false
                break
            end
        end
    end

    for (i, test) in enumerate(tests)
        if test
            push!(A, array_of_sets[i])
        end
    end

    return A
end

function get_one_truefalse_set_from_vector(bool_vector)
    # get the indices of Trues in bool_vector as a BitSet
    # bool vector must be 1-dimensional
    T = Array{Int64,1}()
    F = Array{Int64,1}()

    for i in 1:length(bool_vector)
        if bool_vector[i]
            push!(T, i)
        else
            push!(F, i)
        end
    end

    BT = BitSet(T)
    BF = BitSet(F)

    return BT, BF
end

function get_truefalse_row_indices(bool_array)
    # will return two arrays of BitSets which are:
    # 1) the rows in bool_array that are true, for every column in bool_array - matches between pairs of sequences
    # 2) the rows in bool_array that are false, for every column in bool_array - mismatches between pairs of sequences

    AT = Array{BitSet,1}(undef, size(bool_array, 2))
    AF = Array{BitSet,1}(undef, size(bool_array, 2))

    for i in 1:size(bool_array, 2)
        AT[i], AF[i] = get_one_truefalse_set_from_vector(view(bool_array, :, i))
    end

    return AT, AF
end

function get_sets_in_one_go(the_whole_bool_array)
    # returns an array of bitsets which are the sets

    println("number of threads for set comparisons: ", Threads.nthreads())

    good_views = []
    bad_views = []

    the_whole_thing_is_fine = check_the_view_is_good(the_whole_bool_array)

    if the_whole_thing_is_fine
        return [BitSet(collect(1:size(the_whole_bool_array, 2)))]
    else
        push!(bad_views, the_whole_bool_array)
    end

    while length(bad_views) > 0
        new_bad_views = []
        for bv in bad_views
            new_bitsets = parallel_get_subsets(get_sets_from_view(bv))
            for bitset in new_bitsets
                set = collect(bitset)
                new_view = @view the_whole_bool_array[set, set]
                this_view_is_good = check_the_view_is_good(new_view)
                if this_view_is_good
                    push!(good_views, new_view)
                else
                    push!(new_bad_views, new_view)
                end
            end
        end
        bad_views = new_bad_views
    end

    A = Array{BitSet,1}()

    for gv in good_views
        inds = gv.indices[2]
        push!(A, BitSet(inds))
    end

    return Set(A)
end

# Functions for choosing haplotypes to represent sets
# ---------------------------------------------------

# returns the INDICES of a vector in sorted order. Like sortperm()
# but with ties broken by a second vector
function double_sort_perm(main_vector, secondary_vector;
                          main_increasing = false, secondary_increasing = false)

    @assert length(main_vector) == length(secondary_vector)

    final_order = Array{Int, 1}(undef, length(main_vector))

    sp_1 = sortperm(main_vector, rev = !main_increasing)

    last_item = nothing
    chunk_size = 0
    start = 0
    for (i,item) in enumerate(main_vector[sp_1])
        if i == 1
            last_item = item
            continue
        end

        if item == last_item
            if chunk_size == 0
                start = i - 1
            end

            last_item = item
            chunk_size += 1

        elseif chunk_size > 0
            new_chunk_order = sortperm(secondary_vector[sp_1][start:start + chunk_size], rev = !secondary_increasing)
            sp_1[start:start + chunk_size] = sp_1[start:start + chunk_size][new_chunk_order]
            chunk_size = 0
            last_item = item

        elseif i < length(main_vector[sp_1])
            last_item = item

        end

        if chunk_size > 0
            new_chunk_order = sortperm(secondary_vector[sp_1][start:start + chunk_size], rev = !secondary_increasing)
            sp_1[start:start + chunk_size] = sp_1[start:start + chunk_size][new_chunk_order]
        end

    end

    final_order = sp_1

    return final_order
end

# score sequences based on the number of sequences
# they're different from
function score_differences(pairwise_identity)

    different_scores = Array{Int, 1}(undef, size(pairwise_identity, 2))

    for i in 1:size(pairwise_identity, 2)
        different_scores[i] = size(pairwise_identity, 1) - sum(pairwise_identity[:,i])
    end

    return different_scores
end

# This is the new heuristic - use sequences that are maximally different from all
# other sequences to represent sets. Ties are broken using genome completeness
function write_most_different_unused_seq(pairwise_identity, final_sets, retained, whole_nuc_bit_array, whole_ID_array, alignment_out, append_IDs)

    completeness_scores = score_alignment(whole_nuc_bit_array)
    difference_scores = score_differences(pairwise_identity)

    # convert the set to an Array
    final_sets_as_array = collect(final_sets)

    # split up singletons and others to relieve the burden on sorting?
    singletons = Array{BitSet,1}()
    non_single = Array{BitSet,1}()

    for set in final_sets_as_array
        if length(set) == 1
            push!(singletons, set)
        else
            push!(non_single, set)
        end
    end

    # sort non-singletons by length in increasing order
    sort!(non_single, by = length)

    # and here are singletons + sorted longer sets, concatenated together
    sorted_final_sets = vcat(singletons, non_single)

    # initiate a ditionary that we will use as a check for whether a particular sequence
    # has been used (to represent a set) yet
    used_names = Dict{String, Array{String, 1}}()

    open(alignment_out, "w") do io
        for (i, set) in enumerate(sorted_final_sets)

            # cs is the columns in the alignment that this set represents
            cs = collect(set)

            # If it's a singleton it can't have been so we can just get on with things.
            if length(cs) == 1

                id = ">" * whole_ID_array[cs[1]]
                seq = get_seq_from_1D_byte_array(whole_nuc_bit_array[:,cs[1]])

                # write fasta header and sequence:
                println(io, id)
                println(io, seq)

                used_names[whole_ID_array[cs[1]]] = []

            # otherwise get the most-different-from-all-other-sequences sequence
            # (breaking ties with genome completeness)
            else

                set_completeness_scores = completeness_scores[cs]
                set_difference_scores = difference_scores[cs]
                set_nuc_array = whole_nuc_bit_array[:,cs]
                set_IDs = whole_ID_array[cs]

                max_indices = double_sort_perm(set_difference_scores, set_completeness_scores,
                                               main_increasing = false, secondary_increasing = false)

                indx = 1
                best_column = max_indices[indx]
                best_ID = set_IDs[best_column]

                # A while loop to check that the highest scoring sequence hasn't
                # already been used to represent a set - if it has, go through
                # each next best sequence in turn.
                while in(best_ID, keys(used_names)) && indx < length(max_indices)
                    indx += 1
                    best_column = max_indices[indx]
                    best_ID = set_IDs[best_column]
                end

                # if, after the iteration above, all of the IDs in this set are already
                # representing sets (are in used_names), then we can just skip this set (because
                # all its members are assigned to a set anyway).
                if in(best_ID, keys(used_names))
                    continue
                end

                best_seq = set_nuc_array[:,best_column]
                other_IDs = setdiff(set_IDs, [best_ID])

                used_names[best_ID] = collect(other_IDs)

                if !append_IDs
                    id = ">" * string(best_ID)
                else
                    id = ">" * join(vcat(best_ID, other_IDs), "|")
                end

                seq = get_seq_from_1D_byte_array(best_seq)

                # write fasta header and sequence:
                println(io, id)
                println(io, seq)

            end

        end

    # add retained sequences to the alignment if they aren't included already,
    # and add them to used names so that don't get duplicated in the mapping
    for seqname in retained
        if !in(seqname, keys(used_names))

            indx = findfirst(isequal(seqname), whole_ID_array)

            if indx == nothing
                println(stderr, "warning: " * seqname * " can't be retained because it isn't in the input alignment")
                continue
            end

            if seqname != whole_ID_array[indx]
                throw("bad indexing when forced to retain " * seqname)
            end

            id = ">" * whole_ID_array[indx]
            seq = get_seq_from_1D_byte_array(whole_nuc_bit_array[:,indx])

            println(io, id)
            println(io, seq)

            used_names[seqname] = []
        end
    end

    close(io)
    end

    return used_names
end

# # this is the old heuristic, using the most complete genomes to represent sets
# function write_highest_scoring_unused_seq(final_sets, retained, whole_nuc_bit_array, whole_ID_array, alignment_out, append_IDs)
#
#     all_scores = score_alignment(whole_nuc_bit_array)
#
#     # convert the set to an Array
#     final_sets_as_array = collect(final_sets)
#
#     # split up singletons and others to relieve the burden on sorting?
#     singletons = Array{BitSet,1}()
#     non_single = Array{BitSet,1}()
#
#     for set in final_sets_as_array
#         if length(set) == 1
#             push!(singletons, set)
#         else
#             push!(non_single, set)
#         end
#     end
#
#     # sort non-singletons by length in increasing order
#     sort!(non_single, by = length)
#
#     # and here are singletons + sorted longer sets, concatenated together
#     sorted_final_sets = vcat(singletons, non_single)
#
#     # initiate a ditionary that we will use as a check for whether a particular sequence
#     # has been used (to represent a set) yet
#     used_names = Dict{String, Array{String, 1}}()
#
#     open(alignment_out, "w") do io
#         for (i, set) in enumerate(sorted_final_sets)
#
#             # cs is the columns in the alignment that this set represents
#             cs = collect(set)
#
#             # If it's a singleton it can't have been so we can just get on with things.
#             if length(cs) == 1
#
#                 id = ">" * whole_ID_array[cs[1]]
#                 seq = get_seq_from_1D_byte_array(whole_nuc_bit_array[:,cs[1]])
#
#                 # write fasta header and sequence:
#                 println(io, id)
#                 println(io, seq)
#
#                 used_names[whole_ID_array[cs[1]]] = []
#
#             # otherwise get the highest scoring sequence
#             else
#
#                 set_scores = all_scores[cs]
#                 set_nuc_array = whole_nuc_bit_array[:,cs]
#                 set_IDs = whole_ID_array[cs]
#
#                 max_indices = sortperm(set_scores, rev = true)
#
#                 indx = 1
#                 best_column = max_indices[indx]
#                 best_ID = set_IDs[best_column]
#
#                 # A while loop to check that the highest scoring sequence hasn't
#                 # already been used to represent a set - if it has, go through
#                 # each next best sequence in turn.
#                 while in(best_ID, keys(used_names)) && indx < length(max_indices)
#                     indx += 1
#                     best_column = max_indices[indx]
#                     best_ID = set_IDs[best_column]
#                 end
#
#                 # if, after the iteration above, all of the IDs in this set are already
#                 # representing sets (are in used_names), then we can just skip this set (because
#                 # all it's members are assigned to a set anyway).
#                 if in(best_ID, keys(used_names))
#                     continue
#                 end
#
#                 best_seq = set_nuc_array[:,best_column]
#                 other_IDs = setdiff(set_IDs, [best_ID])
#
#                 used_names[best_ID] = collect(other_IDs)
#
#                 if !append_IDs
#                     id = ">" * string(best_ID)
#                 else
#                     id = ">" * join(vcat(best_ID, other_IDs), "|")
#                 end
#
#                 seq = get_seq_from_1D_byte_array(best_seq)
#
#                 # write fasta header and sequence:
#                 println(io, id)
#                 println(io, seq)
#
#             end
#
#         end
#
#     # add retained sequences to the alignment if they aren't included already,
#     # and add them to used names so that don't get duplicated in the mapping
#     for seqname in retained
#         if !in(seqname, keys(used_names))
#
#             indx = findfirst(isequal(seqname), whole_ID_array)
#
#             if indx == nothing
#                 println(stderr, "warning: " * seqname * " can't be retained because it isn't in the input alignment")
#                 continue
#             end
#
#             if seqname != whole_ID_array[indx]
#                 throw("bad indexing when forced to retain " * seqname)
#             end
#
#             id = ">" * whole_ID_array[indx]
#             seq = get_seq_from_1D_byte_array(whole_nuc_bit_array[:,indx])
#
#             println(io, id)
#             println(io, seq)
#
#             used_names[seqname] = []
#         end
#     end
#
#     close(io)
#     end
#
#     return used_names
# end

function get_redundant_seq_to_tip_relationships(representative_to_set)
    #= representative_to_set is a dict with representative
       haplotypes to be written to file as keys, and the other
       members of the set as values.

       Want to get a map from each possible member => representatives
    =#

    D = Dict{String, Array{String, 1}}()

    for (key, value) in representative_to_set
        # key is the representative sequence for this set
        # value is an array of other sequences in the set
        for v in value
            # if this guy is already a representative of a set,
            # then skip it
            if in(v, keys(representative_to_set))
                continue
            end

            if in(v, keys(D))
                D[v] = vcat(D[v], [key])
            else
                D[v] = [key]
            end
        end
    end

    return D
end

function write_redundant_seq_to_tip_relationships(D, filepath)
    # write a table of sequences that are not included in the
    # collapsed alignment, with information about their place(s)
    # on the tree i.e., for each excluded seq, what tip it can
    # be placed at
    open(filepath, "w") do io
        println(io, "sequence,tips")
        for (key, value) in D
            println(io, key * "," * join(value, "|"))
        end
    close(io)
    end
end

function write_tip_to_redundant_seq_relationships(D, filepath, retained)
    # write the converse file to write_redundant_seq_to_tip_relationships()
    # i.e., tips of the tree and what other seqs they represent (if any)
    open(filepath, "w") do io
        println(io, "tip,sequences")
        for (key, value) in D
            # subtract retained sequences from the set that is represented
            value = setdiff(value, retained)
            # and subtract tips
            value = setdiff(value, keys(D))
            # # ignore sequences that represent no other sequences:
            # if length(value) == 0
            #     continue
            # end
            println(io, key * "," * join(value, "|"))
        end
    close(io)
    end
end

function TEST_get_levels(collection_of_bitsets)
    #=
    return the set of all rows
    that make it into the output - if things
    have worked, then the length of this should
    be the same as the number of sequences - all
    sequences are included in the dataset (at least once)
    =#
    levels = []
    for x in collection_of_bitsets
        for y in x
            push!(levels, y)
        end
    end
    return Set(levels)
end

function TEST_compare_within_sets(array_of_bitsets, nuc_bit_array)
    #=
    compare sequences within sets using the original data -
    =#

    function TEST_all_same(my_small_array)
        for a in 1:size(my_small_array,2) - 1
            for b in (a + 1):size(my_small_array,2)
                for r in 1:size(my_small_array,1)
                    @inbounds x = my_small_array[r,a]
                    @inbounds y = my_small_array[r,b]
                    different = (x & y) < 16

                    # store the different True/False result in an array:
                    if different
                        return false
                    end
                end
            end
        end

        return true
    end

    samples = sample(collect(1:length(array_of_bitsets)), length(array_of_bitsets), replace = false)
    fails = []
    for smp in samples
        bitset = array_of_bitsets[smp]
        columns = collect(bitset)
        if length(columns) > 1
            sub_nuc_bit_array = nuc_bit_array[:,columns]
            test = TEST_all_same(sub_nuc_bit_array)

            if !test
                push!(fails, columns)
            end
        end
    end

    if length(fails) > 0
        return false
    end

    return true
end
