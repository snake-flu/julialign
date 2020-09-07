include("../src/functions_encoding.jl")

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

# this is from functions_encoding.jl:
byte_dict = make_byte_dict()

function test_nucs()
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

    println(all(tests))
end

test_nucs()











#
