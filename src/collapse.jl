using ArgParse

include("functions_align.jl")
include("functions_collapse.jl")
include("functions_collapse_additive.jl")
include("functions_deletions.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--infile", "-i"
            help = "the alignment to collapse, in fasta format"
            required = true
        "--reference"
            help = "the reference sequence, in fasta format"
            required = true
        "--outfile", "-o"
            help = "the name of the reduced alignment to write"
            default = "out.fasta"
        "--rel1"
            help = "output csv with relationships between input sequences, one line per retained sequence"
            default = "tip_to_redundants.csv"
        "--rel2"
            help = "output csv with relationships between input sequences, one line per omitted sequence"
            default = "redundant_to_tips.csv"
        "--retain"
            help = "force the retention of this sequence in the output alignment"
            action = "append_arg"
        "--dels-file"
            help = "file of deletions to type and include as variants when collapsing"
        "--append-as-snp"
            action = "store_true"
            help = "if invoked, and a dels file is provided, append the deletion genotypes as SNPs in the output alignment"
        "--check"
            action = "store_true"
            help = "if invoked, check the sets for not containing any differences"
        "--additive"
            help = "provide a file of previous relationships along with the alignment to collapse
                    - new sequences will be assigned to existing groups if possible, or will form
                    new groups of their own. The full procedure will not be run."
    end

    return parse_args(s)
end

function deletions(parsed_args, ref_array, nucleotide_array, fasta_IDs)
    dels = read_del_file(parsed_args["dels-file"])
    nucleotide_array_dels, A_genotypes_table = type_deletions_and_append_as_SNP(nucleotide_array, dels)

    if Threads.nthreads() == 1
        pairwise_identity = reduced_is_same(nucleotide_array, ref_array)
    else
        pairwise_identity = parallel_reduced_is_same(nucleotide_array, ref_array)
    end

    println("sequence comparisons complete")

    final_sets = get_sets_in_one_go(pairwise_identity)

    levels = TEST_get_levels(final_sets)
    println("number of sequences in input dataset: ", length(fasta_IDs))
    println("number of sequences in output dataset: ", length(levels))

    if parsed_args["check"]
        safe_sets_test = TEST_compare_within_sets(collect(final_sets), nucleotide_array_dels)
    end

    if parsed_args["append-as-snp"]
        tip_to_seq_relationships = get_output_sets(pairwise_identity, final_sets, parsed_args["retain"], nucleotide_array_dels, fasta_IDs)
    else
        tip_to_seq_relationships = get_output_sets(pairwise_identity, final_sets, parsed_args["retain"], nucleotide_array, fasta_IDs)
    end

    write_fasta(tip_to_seq_relationships, nucleotide_array, fasta_IDs, parsed_args["outfile"], false)

    if parsed_args["check"]
        println("safe sets test (true/false): ", safe_sets_test)
    end

    println("number of sequences in output alignment: ", length(keys(tip_to_seq_relationships)))

    write_tip_to_redundant_seq_relationships(tip_to_seq_relationships, parsed_args["rel1"], parsed_args["retain"])

    seq_to_tip_relationships = get_redundant_seq_to_tip_relationships(tip_to_seq_relationships)

    write_redundant_seq_to_tip_relationships(seq_to_tip_relationships, parsed_args["rel2"])
end

function additive(parsed_args, ref_array, nucleotide_array, fasta_IDs)

    old_sets, old_ID_list = read_relationships(parsed_args["additive"])

    old_array, new_array, old_IDs, new_IDs = split_data(nucleotide_array, fasta_IDs, old_ID_list)

    # println("additive mode, number of new sequences to add: ", length(new_IDs))

    if length(new_IDs) == 0
        e = error("no new sequences in alignment - did you mean to use --additive?")
        throw(e)
    end

    m_n_pairwise_identity = m_n_reduced_is_same(new_array, old_array, ref_array)

    new_pairwise_identity = reduced_is_same(new_array, ref_array)

    updated_tip_to_seq_relationships, unplaced_new_ID_indices = assign_new_seqs(m_n_pairwise_identity, new_pairwise_identity, old_IDs, new_IDs, old_sets)

    unplaced_IDs = new_IDs[unplaced_new_ID_indices]

    if length(unplaced_new_ID_indices) > 1

        unplaced_pairwise_identity = new_pairwise_identity[unplaced_new_ID_indices, unplaced_new_ID_indices]

        unplaced_sets = get_sets_in_one_go(unplaced_pairwise_identity)

        unplaced_nuc_bit_array = new_array[:,unplaced_new_ID_indices]

        unplaced_tip_to_seq_relationships = get_output_sets(unplaced_pairwise_identity, unplaced_sets, [], unplaced_nuc_bit_array, unplaced_IDs)

        combined_tip_to_seq_relationships = merge(updated_tip_to_seq_relationships, unplaced_tip_to_seq_relationships)

    elseif length(unplaced_IDs) == 1

        combined_tip_to_seq_relationships = merge(updated_tip_to_seq_relationships, Dict{String, Array{String,1}}(unplaced_IDs[1] => []))

    else

        combined_tip_to_seq_relationships = updated_tip_to_seq_relationships
    end

    write_fasta(combined_tip_to_seq_relationships, nucleotide_array, fasta_IDs, parsed_args["outfile"], false)

    write_tip_to_redundant_seq_relationships(combined_tip_to_seq_relationships, parsed_args["rel1"], [])

    combined_seq_to_tip_relationships = get_redundant_seq_to_tip_relationships(combined_tip_to_seq_relationships)

    write_redundant_seq_to_tip_relationships(combined_seq_to_tip_relationships, parsed_args["rel2"])

end

function vanilla(parsed_args, ref_array, nucleotide_array, fasta_IDs)
    if Threads.nthreads() == 1
        pairwise_identity = reduced_is_same(nucleotide_array, ref_array)
    else
        pairwise_identity = parallel_reduced_is_same(nucleotide_array, ref_array)
    end

    println("sequence comparisons complete")

    final_sets = get_sets_in_one_go(pairwise_identity)

    levels = TEST_get_levels(final_sets)
    println("number of sequences in input dataset: ", length(fasta_IDs))
    println("number of sequences in output dataset: ", length(levels))

    if parsed_args["check"]
        safe_sets_test = TEST_compare_within_sets(collect(final_sets), nucleotide_array)
    end

    tip_to_seq_relationships = get_output_sets(pairwise_identity, final_sets, parsed_args["retain"], nucleotide_array, fasta_IDs)

    write_fasta(tip_to_seq_relationships, nucleotide_array, fasta_IDs, parsed_args["outfile"], false)

    if parsed_args["check"]
        println("safe sets test (true/false): ", safe_sets_test)
    end

    println("number of sequences in output alignment: ", length(keys(tip_to_seq_relationships)))

    write_tip_to_redundant_seq_relationships(tip_to_seq_relationships, parsed_args["rel1"], parsed_args["retain"])

    seq_to_tip_relationships = get_redundant_seq_to_tip_relationships(tip_to_seq_relationships)

    write_redundant_seq_to_tip_relationships(seq_to_tip_relationships, parsed_args["rel2"])
end

function main()
    parsed_args = parse_commandline()

    ref_array, ref_id = populate_byte_array_get_names(parsed_args["reference"])
    nucleotide_array, fasta_IDs = populate_byte_array_get_names(parsed_args["infile"])

    if length(ref_array) != size(nucleotide_array, 1)
        e = error("reference and alignment do not have the same number of sites")
        throw(e)
    end

    if parsed_args["additive"] != nothing
        if parsed_args["dels-file"] != nothing
            e = error("can't use a dels file with --additive at the moment")
            throw(e)
        end

        additive(parsed_args, ref_array, nucleotide_array, fasta_IDs)

    elseif parsed_args["dels-file"] != nothing

        deletions(parsed_args, ref_array, nucleotide_array, fasta_IDs)

    else

        vanilla(parsed_args, ref_array, nucleotide_array, fasta_IDs)

    end
end

main()
