# byte_dict is a map from the Char representation of IUPAC nucleotide
# codes to uint8 values according to the (modified) bit-level coding scheme
# for nucleotides developed by Emmanuel Paradis:
# http://ape-package.ird.fr/misc/BitLevelCodingScheme.html
function make_byte_dict()
    byte_dict = Dict{Char, UInt8}('A' => 136,
                                    'G' => 72,
                                    'C' => 40,
                                    'T' => 24,
                                    'R' => 192,
                                    'M' => 160,
                                    'W' => 144,
                                    'S' => 96,
                                    'K' => 80,
                                    'Y' => 48,
                                    'V' => 224,
                                    'H' => 176,
                                    'D' => 208,
                                    'B' => 112,
                                    'N' => 240,
                                    '-' => 244,
                                    '?' => 242)

    return byte_dict
end


# nuc_dict maps from uint8 representation of nucleotides to their IUPAC code
function make_nuc_dict()
    nuc_dict = Dict{UInt8, Char}(136 => 'A',
                                 72 => 'G',
                                 40 => 'C',
                                 24 => 'T',
                                 192 => 'R',
                                 160 => 'M',
                                 144 => 'W',
                                 96 => 'S',
                                 80 => 'K',
                                 48 => 'Y',
                                 224 => 'V',
                                 176 => 'H',
                                 208 => 'D',
                                 112 => 'B',
                                 240 => 'N',
                                 244 => '-',
                                 242 => '?')

    return nuc_dict
end


# score_dict is a map from uint8 values for IUPAC nucleotide codes to an
# integer for how unambiguous they are. The score is caculated as:
# 12 * 1/possible real nucleotides. E.g. an 'A' / 136::uint8 scores 12, but an
# 'N' (A, C, G or T) / 240::uint8 scores 3
function make_score_dict()
    score_dict = Dict{UInt8, Int64}(136 => 12,
                                     72 => 12,
                                     40 => 12,
                                     24 => 12,
                                     192 => 6,
                                     160 => 6,
                                     144 => 6,
                                     96 => 6,
                                     80 => 6,
                                     48 => 6,
                                     224 => 4,
                                     176 => 4,
                                     208 => 4,
                                     112 => 4,
                                     240 => 3,
                                     244 => 3,
                                     242 => 3)

    return score_dict
end
