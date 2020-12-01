using ArgParse

include("functions_closest.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--target"
            help = "The alignment to search for closest matches, in fasta format"
            required = true
        "--query"
            help = "The sequence(s) to find matches for, in fasta format"
            required = true
        "--reference", "-r"
            help = "The reference sequence, in fasta format"
            required = false
        "--csv-out", "-o"
            help = "The name of csv file containing closest matches to write"
            default = "out.csv"
        "--fast"
            help = "Restrict comparison to sites that are different from the reference.
                    Requires --reference to be provided; doesn't calculate a per-site distance"
            action = :store_true
    end

    return parse_args(s)
end

function main()

    parsed_args = parse_commandline()

    if parsed_args["fast"]
        if parsed_args["reference"] == nothing
            e = error("need to provide --reference if you use --fast")
            throw(e)
        end

        fast_closest(parsed_args["target"], parsed_args["query"], parsed_args["reference"], parsed_args["csv-out"])
    else
        closest(parsed_args["target"], parsed_args["query"], parsed_args["csv-out"])
    end
end

main()
