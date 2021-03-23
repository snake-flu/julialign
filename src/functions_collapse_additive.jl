# include("functions_align.jl")
# include("functions_collapse.jl")

### --------- functions for additive addition of new seqs

# struct redundancy_set
#     tip::String
#     redundants::Array{String, 1}
# end

function m_n_reduced_is_same(nuc_bit_array_A, nuc_bit_array_B, ref_array)

    width = size(nuc_bit_array_A, 2)
    height = size(nuc_bit_array_B, 2)

    same_array = trues(height, width)

    # println("number of threads for sequence comparisons: ", Threads.nthreads())

    difference_list_A = get_list_of_differences(nuc_bit_array_A, ref_array)
    difference_list_B = get_list_of_differences(nuc_bit_array_B, ref_array)

    for a in 1:size(nuc_bit_array_A, 2)

        for b in 1:size(nuc_bit_array_B, 2)

            start_again = true

            for r in difference_list_A[a]
                @inbounds x = nuc_bit_array_A[r,a]
                @inbounds y = nuc_bit_array_B[r,b]
                different = (x & y) < 16

                # store the different True/False result in an array:
                if different
                    same_array[b,a] = !different
                    start_again = false
                    break
                end
            end

            # need to check the second seq's sites too, if we get this far:
            if start_again
                for r in difference_list_B[b]
                    if r in difference_list_A[a]
                        continue
                    end

                    @inbounds x = nuc_bit_array_A[r,a]
                    @inbounds y = nuc_bit_array_B[r,b]
                    different = (x & y) < 16

                    # store the different True/False result in an array:
                    if different
                        same_array[b,a] = !different
                        break
                    end
                end
            end
        end
    end

    return same_array

end

function read_relationships(filepath)
    # filepath should point to a tip_2_redundants.csv format file

    all_seqs = Array{String,1}()
    sets = Dict{String, Array{String, 1}}()

    first = true

    open(filepath, "r") do io
        while !eof(io)

            if first
                readline(io)
                first = false
                continue
            end

            line = readline(io)

            fields = split(line, ",")
            tip = fields[1]

            push!(all_seqs, tip)

            if length(fields) != 2
                e = error("badly formed additions file")
                throw(e)
            elseif fields[2] == ""
                redundants = Array{String,1}()
            else
                redundants = split(fields[2], "|")
                for redundant in redundants
                    push!(all_seqs, redundant)
                end
            end

            sets[tip] = redundants
        end
    end

    all_seqs = Set(all_seqs)

    return sets, all_seqs
end

function update_relationships_based_on_appearance(sets, all_seqs, fasta_IDs)
    #=
    sets = a dict of string(tip):array of strings(redundants)
    all_seqs = a set of strings (all the sequences from the relationships file)
    fasta_IDs = an array of strings which are the IDs for all the seqs in the fasta file

    A function to update the relationships that have been read in from the
    previous "tip_to_redundants.csv" file, to deal with the case where sequences
    that were in the previous dataset are missing in the new dataset.

    The following things should happen:

    1) for sets with tips that are not in the new dataset, a new tip should be chosen
       from within that set. If this is the case this redundant also needs to be removed
       from any other sets where it is a redundant (as would happen in the full routine).
       If the set has no other sequences, it can be thrown out.
    2) for sets with redundants that are not in the new dataset, those redundants
       should be stripped out of the sets.
    =#

    # make a set of the fasta IDs for faster lookup
    fasta_ID_set = Set(fasta_IDs)

    # sets whose tips is missing in these data. will check that all their
    # redundants are not redundants in sets_that_dont_need_new_tips already
    # (unchanged from the input data), because if so, can throw out this
    # whole set
    sets_that_need_new_tips = Dict{String, Array{String, 1}}()
    sets_that_dont_need_new_tips = Dict{String, Array{String, 1}}()

    # make a set of redundants that have become tips and therefore need to be taken
    # out of any other sets that they are in.

    # step 1: sort out sets that need new tips (or could potentially be thrown away)
    for (key, value) in sets
        # if this key doesn't exist any more
        if !(key in fasta_ID_set)
            # there are any redundants, flag for followup (else ignore)
            if length(value) >= 1
                sets_that_need_new_tips[key] = value
            end
        else
            sets_that_dont_need_new_tips[key] = value
        end
    end

    # step 2:
    # a) get the set of all the redundant seqs in the unchanged data:
    # make a set of the current redundants for faster lookup
    unchanged_set_seqs = Set{String}()
    for (key, value) in sets_that_dont_need_new_tips
        for v in value
            if v in fasta_ID_set
                push!(unchanged_set_seqs, v)
            end
        end
    end

    # b) find new tips for sets whose tip is missing, or throw them out if
    # they're already represented
    redundants_that_are_now_tips = Set{String}()
    sets_with_new_tips = Dict{String, Array{String, 1}}()
    for (key, value) in sets_that_need_new_tips
        # key doesn't exist in the new dataset
        for v in value
            # v might not:
            if !(v in fasta_ID_set)
                continue
            # or it might already be represented elsewhere:
            elseif v in unchanged_set_seqs
                continue
            # otherwise just
            else
                sets_with_new_tips[v] = setdiff(value, v)
                push!(redundants_that_are_now_tips, v)
                break
            end
        end
    end

    # step 3: merge the unchanged and the new sets
    new_sets_temp = Dict{String, Array{String, 1}}()

    for (key, value) in sets_with_new_tips
        new_sets_temp[key] = value
    end

    for (key, value) in sets_that_dont_need_new_tips
        new_sets_temp[key] = value
    end

    # step 4: remove redundants that are missing or have become tips
    # and populate the new set of old seqs at the same time. ???

    new_sets = Dict{String, Array{String, 1}}()
    new_all_seqs = Set{String}()

    for (key, value) in new_sets_temp
        push!(new_all_seqs, key)

        newredundants = Array{String, 1}()
        for redundant in value
            if !(redundant in redundants_that_are_now_tips)
                if redundant in fasta_ID_set
                    push!(newredundants, redundant)
                    push!(new_all_seqs, redundant)
                end
            end
        end

        new_sets[key] = newredundants
    end

    return new_sets, new_all_seqs
end

function update_relationships_based_on_sequences(sets, all_seqs, nucleotide_array, fasta_IDs, ref_array)
    #=
    sets = a dict of string(tip):array of strings(redundants)
    all_seqs = a set of strings (all the sequences from the relationships file)
    nucleotide_array = the alignment
    fasta_IDs = an array of strings which are the IDs for all the seqs in the fasta file

    A function to update the relationships that have been read in from the
    previous "tip_to_redundants.csv" file, to deal with the case where sequences
    that were in the previous dataset are different in the new dataset (the only
    way you can tell is by testing the previous sets for still being good sets).

    check that the sets are still good sets. TODO im not sure how to deal with the situation
    where they are not. bad apples (sequences) can be pushed into the new set, which means that
    they will get reassigned additively, but im not sure how to identify the bad apples in
    the first place. I could throw the whole set out. - this might be easiest?
    =#

    # need to subset the nucleotide array, so need to know what columns to use:
    ID_to_column_dict = get_id_to_column_dict(fasta_IDs) # returns a dict fasta IDs to their column in the alignment

    new_sets = Dict{String, Array{String, 1}}()
    new_all_seqs = Set{String}()

    for (key, value) in sets
        # an array containing all the seq names for this set:
        IDs = Array{String, 1}()
        push!(IDs, key)
        for v in value
            push!(IDs, v)
        end

        # an array containing all the column numbers for this set:
        colmns = Array{Int32, 1}()
        for ID in IDs
            push!(colmns, ID_to_column_dict[ID])
        end

        # now check the set:
        same_array = reduced_is_same(view(nucleotide_array, :, colmns), ref_array) # returns the same_array

        set_is_good = check_the_view_is_good(same_array) # returns true/false the same_array is ok

        if set_is_good
            # if the set is still fine, we can keep it
            new_sets[key] = value

            # and we add its members to the set of all sequences
            push!(new_all_seqs, key)
            for v in value
                push!(new_all_seqs, v)
            end
        else
            # Otherwise, we don't keep it. But we don't need to strip its members from other sets,
            # because they will be tested anyway...
            continue
        end

    end

    return new_sets, new_all_seqs
end

#
function split_data(nuc_bit_array, fasta_ID_array, old_IDs)

    old_indices = Array{Int, 1}()
    sizehint!(old_indices, length(old_IDs))
    new_indices = Array{Int, 1}()
    sizehint!(new_indices, length(old_IDs))

    for (i, ID) in enumerate(fasta_ID_array)
        if ID in old_IDs
            push!(old_indices, i)
        else
            push!(new_indices, i)
        end
    end

    old_array = nuc_bit_array[:,old_indices]
    new_array = nuc_bit_array[:,new_indices]

    old_IDs = fasta_ID_array[old_indices]
    new_IDs = fasta_ID_array[new_indices]

    return old_array, new_array, old_IDs, new_IDs
end

# function get_set_dict(old_sets)
#     set_dict = Dict{String, Array{String, 1}}()
#     for set in old_sets
#         set_dict[set.tip] = set.redundants
#     end
#
#     return set_dict
# end

# # comparison between two vectors (at sites that are diff from ref in either)
# function check_identity(nuc_bit_V_A, nuc_bit_V_B, diff_list_A, diff_list_B)
#
#     start_again = true
#     all_same = true
#
#     for r in diff_list_A
#         @inbounds x = nuc_bit_V_A[r]
#         @inbounds y = nuc_bit_V_B[r]
#         different = (x & y) < 16
#
#         # store the different True/False result in an array:
#         if different
#             all_same = false
#             start_again = false
#             break
#         end
#     end
#
#     # need to check the second seq's sites too, if we get this far:
#     if start_again
#         for r in diff_list_B
#             if r in diff_list_A
#                 continue
#             end
#
#             @inbounds x = nuc_bit_V_A[r]
#             @inbounds y = nuc_bit_V_B[r]
#             different = (x & y) < 16
#
#             # store the different True/False result in an array:
#             if different
#                 all_same = false
#                 break
#             end
#         end
#     end
#
#     return all_same
# end

function assign_new_seqs(m_n_pairwise_identity, unplaced_pairwise_identity, old_IDs, new_IDs, old_sets)

    old_ID_set = Set(old_IDs)

    # old_set_dict = get_set_dict(old_sets)
    old_ID_to_old_cl_dict = get_id_to_column_dict(old_IDs)
    new_ID_to_new_cl_dict = get_id_to_column_dict(new_IDs)

    placed_new_tips = Array{String, 1}()

    for (i, ID) in enumerate(new_IDs)

        for tip in keys(old_sets)

            tip_column = old_ID_to_old_cl_dict[tip]

            tip_same = m_n_pairwise_identity[tip_column, i]

            if tip_same
                 # check all the other members of the set...

                 set_same = true

                 other_members = old_sets[tip]

                 for member in other_members
                     # possible that this member is a new one (that has already been added to this set),
                     # so need to check what to do here:
                     if member in old_ID_set
                         redundant_column = old_ID_to_old_cl_dict[member]
                         set_same = m_n_pairwise_identity[redundant_column, i]
                         if ! set_same
                             break
                         end
                     else
                         # CHANGE BELOW - TO LOOK UP IN NEW NUC ARRAY
                         redundant_column = new_ID_to_new_cl_dict[member]
                         set_same = unplaced_pairwise_identity[redundant_column, i]
                         if ! set_same
                             break
                         end
                     end
                 end

                 if set_same
                     # update the sets
                     old_sets[tip] = vcat(old_sets[tip], ID)
                     #
                     push!(placed_new_tips, ID)
                     #
                 end
            end

        end
    end

    placed_new_tips = Set(placed_new_tips)

    unplaced_new_tip_indices = Array{Int, 1}()
    for (i, tip) in enumerate(new_IDs)
        if ! (tip in placed_new_tips)
            push!(unplaced_new_tip_indices, i)
        end
    end

    updated_sets = old_sets

    return updated_sets, unplaced_new_tip_indices
end

function TEST_get_dict_levels(relationships_dict)
    s = Set{String}()
    for (key, value) in relationships_dict
        push!(s, key)
        for v in value
            push!(s, v)
        end
    end
    return s
end

function TEST_get_sets_from_relationships(fasta_IDs, relationships_dict)

    col_dict = get_id_to_column_dict(fasta_IDs)

    A = Array{BitSet,1}()

    for (key, value) in relationships_dict
        temp_A = BitSet()
        push!(temp_A, col_dict[key])
        for v in value
            push!(temp_A, col_dict[v])
        end
        
        push!(A, temp_A)
    end

    return A
end

#
# function assign_new_seqs(old_nuc_bit_array, new_nuc_bit_array, old_IDs, new_IDs, old_sets, ref_array)
#
#     # old_set_dict = get_set_dict(old_sets)
#     old_ID_to_old_cl_dict = get_id_to_column_dict(old_IDs)
#     new_ID_to_new_cl_dict = get_id_to_column_dict(new_IDs)
#
#     difference_list_new = get_list_of_differences(new_nuc_bit_array, ref_array)
#     difference_list_old = get_list_of_differences(old_nuc_bit_array, ref_array)
#
#     # updated_sets = Array{redundancy_set, 1}()
#
#     placed_new_tips = Array{String, 1}()
#
#     for (i, ID) in enumerate(new_IDs)
#
#         for tip in keys(old_sets)
#
#             tip_column = old_ID_to_old_cl_dict[tip]
#
#             tip_same = check_identity(old_nuc_bit_array[:,tip_column], new_nuc_bit_array[:,i], difference_list_old[tip_column], difference_list_new[i])
#
#             if tip_same
#                  # check all the other members of the set...
#
#                  set_same = true
#
#                  other_members = old_sets[tip]
#
#                  for member in other_members
#                      # possible that this member is a new one (that has already been added to this set),
#                      # so need to check what to do here:
#                      if member in old_IDs
#                          redundant_column = old_ID_to_old_cl_dict[member]
#                          set_same = check_identity(old_nuc_bit_array[:,redundant_column], new_nuc_bit_array[:,i], difference_list_old[redundant_column], difference_list_new[i])
#                          if ! set_same
#                              break
#                          end
#                      else
#                          # CHANGE BELOW - TO LOOK UP IN NEW NUC ARRAY
#                          redundant_column = new_ID_to_new_cl_dict[member]
#                          set_same = check_identity(new_nuc_bit_array[:,redundant_column], new_nuc_bit_array[:,i], difference_list_new[redundant_column], difference_list_new[i])
#                          if ! set_same
#                              break
#                          end
#                      end
#                  end
#
#                  if set_same
#                      # update the sets
#                      old_sets[tip] = vcat(old_sets[tip], ID)
#                      #
#                      push!(placed_new_tips, ID)
#                      #
#                  end
#             end
#
#         end
#     end
#
#     placed_new_tips = Set(placed_new_tips)
#
#     unplaced_new_tip_indices = Array{Int, 1}()
#     for (i,tip) in enumerate(new_IDs)
#         if ! (tip in placed_new_tips)
#             push!(unplaced_new_tip_indices, i)
#         end
#     end
#
#     updated_sets = old_sets
#
#     unplaced_nuc_bit_array = new_nuc_bit_array[:,unplaced_new_tip_indices]
#     unplaced_IDs = new_IDs[unplaced_new_tip_indices]
#
#     return updated_sets, unplaced_nuc_bit_array, unplaced_IDs
# end
