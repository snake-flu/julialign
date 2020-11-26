include("../src/functions_align.jl")
include("../src/functions_closest.jl")
include("../src/functions_collapse.jl")
include("../src/functions_collapse_additive.jl")
include("../src/functions_deletions.jl")
include("../src/functions_encoding.jl")

#=
----------------- encoding tests --------------------
=#

lookup_char = Dict{Char, Array{Char,1}}('A' => ['A'],
                                        'C' => ['C'],
                                        'G' => ['G'],
                                        'T' => ['T'],
                                        'R' => ['A', 'G'],
                                        'Y' => ['C', 'T'],
                                        'S' => ['G', 'C'],
                                        'W' => ['A', 'T'],
                                        'K' => ['G', 'T'],
                                        'M' => ['A', 'C'],
                                        'B' => ['C', 'G', 'T'],
                                        'D' => ['A', 'G', 'T'],
                                        'H' => ['A', 'C', 'T'],
                                        'V' => ['A', 'C', 'G'],
                                        'N' => ['A', 'C', 'G', 'T'],
                                        '?' => ['A', 'C', 'G', 'T'],
                                        '-' => ['A', 'C', 'G', 'T']  )

# this function tests the encoding from char -> paradis uint8
function test_nucs()
    # this is from functions_encoding.jl:
    byte_dict = make_byte_dict()

    nucs = ['A', 'G', 'C', 'T', 'R', 'M', 'W', 'S', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '-', '?']

    tests = []

    for i in eachindex(nucs)
        for j in eachindex(nucs)

            nuc1 = nucs[i]
            nuc2 = nucs[j]

            nuc1_chars = lookup_char[nuc1]
            nuc2_chars = lookup_char[nuc2]

            byte1 = byte_dict[nuc1]
            byte2 = byte_dict[nuc2]

            byte_different = (byte1 & byte2) < 16
            byte_same = (byte1 & 8 == 8) && byte1 == byte2

            nuc_different = length(intersect(nuc1_chars, nuc2_chars)) == 0
            nuc_same = in(nuc1, ['A', 'C', 'G', 'T']) && nuc1 == nuc2

            push!(tests, byte_different == nuc_different && byte_same == nuc_same)
        end
    end

    println("test_nucs() passes test: " * string(all(tests)))
end

test_nucs()

# this function tests the encoding from ascii byte -> paradis uint8
function test_nucs2()
    # this is from functions_encoding.jl:
    byte_array = make_byte_array()

    nucs = ['A', 'G', 'C', 'T', 'R', 'M', 'W', 'S', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '-', '?']

    tests = []

    for i in eachindex(nucs)
        for j in eachindex(nucs)

            nuc1 = nucs[i]
            nuc2 = nucs[j]

            nuc1_chars = lookup_char[nuc1]
            nuc2_chars = lookup_char[nuc2]

            byte1 = byte_array[Int(nuc1)]
            byte2 = byte_array[Int(nuc2)]

            byte_different = (byte1 & byte2) < 16
            byte_same = (byte1 & 8 == 8) && byte1 == byte2

            nuc_different = length(intersect(nuc1_chars, nuc2_chars)) == 0
            nuc_same = in(nuc1, ['A', 'C', 'G', 'T']) && nuc1 == nuc2

            push!(tests, byte_different == nuc_different && byte_same == nuc_same)
        end
    end

    println("test_nucs2() passes test: " * string(all(tests)))
end

test_nucs2()

#=
----------------- fastaio / encodingdeencoding test --------------------
=#

function test_fastaread()

    tests = Vector{Bool}()

    array, ids = populate_byte_array_get_names(dirname(Base.source_path()) * "/read_test.fasta")

    # test the fasta headers are parsed correctly
    id_test = ids == String["A/1", "A/2", "B/1", "C/1", "C/2", "C/3"]
    push!(tests, id_test)

    # test the sequences are encoded and deencoded correctly
    seqs = Array{String, 1}(undef, size(array, 2))
    for i in 1:size(array, 2)
        seqs[i] = get_seq_from_1D_byte_array(array[:,i])
    end

    seq_test = seqs == ["ACGTACGTACGTACGTACGT",
                        "ACGTACGTACGTACGTANNN",
                        "AAGTACGTACGTACGTACGT",
                        "ACCTACGTACGTACGTACGT",
                        "ACCTACGTACGTACGTACNN",
                        "ACGTRMWSKYVHDBN-?ATG"]

    push!(tests, seq_test)

    println("test_fastaio() passes test: " * string(all(tests)))
end

test_fastaread()


#=
----------------- collapse test --------------------
=#

function test_collapse()

    tests = Vector{Bool}()

    nucleotide_array, fasta_IDs = populate_byte_array_get_names(dirname(Base.source_path()) * "/collapse_test.fasta")
    ref_array, ref_ID = populate_byte_array_get_names(dirname(Base.source_path()) * "/collapse_test_ref.fasta")
    pairwise_identity = reduced_is_same(nucleotide_array, ref_array)
    final_sets = get_sets_in_one_go(pairwise_identity)

    safe_sets_test = TEST_compare_within_sets(collect(final_sets), nucleotide_array)
    push!(tests, safe_sets_test)

    tip_to_seq_relationships = get_output_sets(pairwise_identity, final_sets, [], nucleotide_array, fasta_IDs)
    seq_to_tip_relationships = get_redundant_seq_to_tip_relationships(tip_to_seq_relationships)

    tip_to_seq_relationships_known = Dict{String, Array{String, 1}}()
    tip_to_seq_relationships_known["seq8"] = ["seq10", "seq13"]
    tip_to_seq_relationships_known["seq1"] = ["seq5", "seq6", "seq13", "seq14"]
    tip_to_seq_relationships_known["seq4"] = ["seq9", "seq10", "seq13"]
    tip_to_seq_relationships_known["seq3"] = ["seq5", "seq6", "seq13"]
    tip_to_seq_relationships_known["seq12"] = ["seq13"]
    tip_to_seq_relationships_known["seq2"] = ["seq5", "seq7", "seq13"]
    tip_to_seq_relationships_known["seq11"] = ["seq13"]

    tip_to_seq_key_test = keys(tip_to_seq_relationships) == keys(tip_to_seq_relationships_known)
    push!(tests, tip_to_seq_key_test)
    if tip_to_seq_key_test
        for key in keys(tip_to_seq_relationships_known)
            subtest = Set(tip_to_seq_relationships[key]) == Set(tip_to_seq_relationships_known[key])
            push!(tests, subtest)
        end
    end

    seq_to_tip_relationships_known = Dict{String, Array{String, 1}}()
    seq_to_tip_relationships_known["seq14"] = ["seq1"]
    seq_to_tip_relationships_known["seq6"] = ["seq1", "seq3"]
    seq_to_tip_relationships_known["seq5"] = ["seq1", "seq3", "seq2"]
    seq_to_tip_relationships_known["seq9"] = ["seq4"]
    seq_to_tip_relationships_known["seq13"] = ["seq8", "seq1", "seq4", "seq3", "seq12", "seq2", "seq11"]
    seq_to_tip_relationships_known["seq10"] = ["seq8", "seq4"]
    seq_to_tip_relationships_known["seq7"] = ["seq2"]

    seq_to_tip_key_test = keys(seq_to_tip_relationships) == keys(seq_to_tip_relationships_known)
    push!(tests, seq_to_tip_key_test)
    if seq_to_tip_key_test
        for key in keys(seq_to_tip_relationships_known)
            subtest = Set(seq_to_tip_relationships[key]) == Set(seq_to_tip_relationships_known[key])
            push!(tests, subtest)
        end
    end

    println("test_collapse() passes test: " * string(all(tests)))
end

test_collapse()


#=
----------------- score_alignment test --------------------
=#

function test_score_alignment()
    tests = Vector{Bool}()

    array, ids = populate_byte_array_get_names(dirname(Base.source_path()) * "/read_test.fasta")

    known_scores = [240, 213, 240, 240, 222, 145]

    # score1 = score_alignment(array)
    score2 = score_alignment2(array)

    # push!(tests, score1 == known_scores)
    push!(tests, score2 == known_scores)

    println("test_score_alignment() passes test: " * string(all(tests)))
end

test_score_alignment()

#=
----------------- closest test --------------------
=#

function test_closest()

    tests = Vector{Bool}()

    R_A, R_names = populate_byte_array_get_names(dirname(Base.source_path()) * "/closest_test_ref.fasta")
    T_A, T_names = populate_byte_array_get_names(dirname(Base.source_path()) * "/closest_test_target.fasta")
    Q_A, Q_names = populate_byte_array_get_names(dirname(Base.source_path()) * "/closest_test_query.fasta")

    target_completeness = score_alignment2(T_A)

    test_target_completeness = target_completeness == [36, 30, 36, 30, 27]
    push!(tests, test_target_completeness)

    distance = get_difference_matrix(T_A, Q_A)
    test_distance = distance == reshape([0.0/3.0, 0.0/2.0, 1.0/3.0, 1.0/2.0, 0.0/2.0,
                                         0.0/3.0, 0.0/2.0, 1.0/3.0, 1.0/2.0, 0.0/2.0,
                                         1.0/1.0, 1.0/1.0, 1.0/1.0, 1.0/1.0, 1.0/1.0], 5, 3)
    push!(tests, test_distance)

    SNP_distance = get_SNP_distance_matrix(T_A, Q_A, R_A)
    test_SNP_distance = SNP_distance == reshape([0, 0, 1, 1, 0,
                                                 0, 0, 1, 1, 0,
                                                 1, 1, 1, 1, 1], 5, 3)
    push!(tests, test_SNP_distance)

    best_indices_distance = [1, 1, 1]
    best_targets_distance = ["target1", "target1", "target1"]
    best_SNPs_distance = ["", "", "3AG"]

    for query_index in 1:size(distance, 2)
        best_indx = get_best_target_index(distance[:,query_index], target_completeness)
        push!(tests, best_indx == best_indices_distance[query_index])

        SNPs = get_SNPs(Q_A[:,query_index], T_A[:,best_indx])
        push!(tests, SNPs == best_SNPs_distance[query_index])

        Q_name = Q_names[query_index]

        T_name = T_names[best_indx]
        push!(tests, T_name == best_targets_distance[best_indx])
    end

    best_indices_SNPs = [1, 1, 1]
    best_targets_SNPs = ["target1", "target1", "target1"]
    best_SNPs_SNPs = ["", "", "3AG"]

    for query_index in 1:size(distance, 2)
        best_indx = get_best_target_index(SNP_distance[:,query_index], target_completeness)
        push!(tests, best_indx == best_indices_SNPs[query_index])

        SNPs = get_SNPs(Q_A[:,query_index], T_A[:,best_indx])
        push!(tests, SNPs == best_SNPs_SNPs[query_index])

        Q_name = Q_names[query_index]

        T_name = T_names[best_indx]
        push!(tests, T_name == best_targets_SNPs[best_indx])
    end

    println("test_closest() passes test: " * string(all(tests)))
end

test_closest()









#
