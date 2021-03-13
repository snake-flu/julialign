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

    @time ref_array, ref_ID = populate_byte_array_get_names(parsed_args["reference"])

    @time nucleotide_array_1, fasta_IDs = populate_byte_array_get_names(parsed_args["infile1"])

    @time snpmat = pairsnp(nucleotide_array_1, ref_array)

    @time write_pairsnp(parsed_args["outfile"], snpmat, fasta_IDs)

end

main()
