using ArgParse

include("functions_align.jl")
include("functions_pairsnp.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--infile1"
            help = "alignment 1"
            required = true
        "--infile2"
            help = "alignment 2"
        "--reference", "-r"
            help = "the reference sequence"
            required = true
        "--outfile", "-o"
            help = "the name of the file to write"
            default = "pairsnp.csv"
        end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    ref_array, ref_ID = populate_byte_array_get_names(parsed_args["reference"])

    nucleotide_array_1, fasta_IDs_1 = populate_byte_array_get_names(parsed_args["infile1"])

    if parsed_args["infile2"] != nothing
        nucleotide_array_2, fasta_IDs_2 = populate_byte_array_get_names(parsed_args["infile2"])
        snpmat = pairsnp_rectangle(nucleotide_array_1, nucleotide_array_2, ref_array)
        write_pairsnp_rectangle(parsed_args["outfile"], snpmat, fasta_IDs_1, fasta_IDs_2)
    else
        snpmat = pairsnp_square(nucleotide_array_1, ref_array)
        write_pairsnp_square(parsed_args["outfile"], snpmat, fasta_IDs_1)
    end
end

main()
