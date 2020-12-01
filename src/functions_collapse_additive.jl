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
    # filepath should point to a tips_2_redundants.csv format file

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

    return(sets, all_seqs)
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
