using ArgParse

include("functions_align.jl")
include("functions_collapse.jl")
include("functions_deletions.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--infile", "-i"
            help = "the alignment to collapse, in fasta format"
            required = true
        "--reference"
            help = "the reference sequence, in fasta format (for speed?)"
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
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    fasta_IDs = get_fasta_descriptions(parsed_args["infile"])
    ref_array = populate_byte_array(parsed_args["reference"])
    nucleotide_array = populate_byte_array(parsed_args["infile"])

    if parsed_args["dels-file"] != nothing
        dels = read_del_file(parsed_args["dels-file"])
        nucleotide_array_dels, A_genotypes_table = type_deletions_and_append_as_SNP(nucleotide_array, dels)

        if Threads.nthreads() == 1
            # pairwise_identity = is_same(nucleotide_array)
            pairwise_identity = reduced_is_same(nucleotide_array, ref_array)
        else
            # pairwise_identity = parallel_is_same(nucleotide_array_dels)
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
            tip_to_seq_relationships = write_most_different_unused_seq(pairwise_identity, final_sets, parsed_args["retain"], nucleotide_array_dels, fasta_IDs, parsed_args["outfile"], false)
        else
            tip_to_seq_relationships = write_most_different_unused_seq(pairwise_identity, final_sets, parsed_args["retain"], nucleotide_array, fasta_IDs, parsed_args["outfile"], false)
        end

    else
        if Threads.nthreads() == 1
            # pairwise_identity = is_same(nucleotide_array)
            pairwise_identity = reduced_is_same(nucleotide_array, ref_array)
        else
            # pairwise_identity = parallel_is_same(nucleotide_array_dels)
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

        tip_to_seq_relationships = write_most_different_unused_seq(pairwise_identity, final_sets, parsed_args["retain"], nucleotide_array, fasta_IDs, parsed_args["outfile"], false)
    end

    if parsed_args["check"]
        println("safe sets test (true/false): ", safe_sets_test)
    end
    
    println("number of sequences in output alignment: ", length(keys(tip_to_seq_relationships)))

    write_tip_to_redundant_seq_relationships(tip_to_seq_relationships, parsed_args["rel1"], parsed_args["retain"])

    seq_to_tip_relationships = get_redundant_seq_to_tip_relationships(tip_to_seq_relationships)
    write_redundant_seq_to_tip_relationships(seq_to_tip_relationships, parsed_args["rel2"])
end

main()
