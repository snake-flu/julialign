using ArgParse

include("functions_align.jl")
include("functions_fourgametetest.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--infile", "-i"
            help = "alignment"
            required = true
        end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    nucleotide_array, fasta_IDs = populate_byte_array_get_names(parsed_args["infile"])
    get_biallelic_columns(nucleotide_array, fasta_IDs)
end

main()
