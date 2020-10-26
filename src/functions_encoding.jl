#=
See - http://ape-package.ird.fr/misc/BitLevelCodingScheme.html for the bit-level
coding scheme :

Nucleotide	    IUPAC code	 Bit-level code	  Value
---------------------------------------------------
A	                  A	      10001000	       136
G	                  G	      01001000	        72
C	                  C	      00101000	        40
T	                  T	      00011000	        24
A or G	              R	      11000000	       192
A or C	              M	      10100000	       160
A or T	              W	      10010000	       144
G or C	              S	      01100000	        96
G or T	              K	      01010000	        80
C or T	              Y	      00110000	        48
A or G or C	          V	      11100000	       224
A or C or T	          H	      10110000	       176
A or G or T	          D	      11010000	       208
G or C or T	          B	      01110000	       112
A or G or C or T      N	      11110000	       240
Alignment gap        (–)	  00000100	         4
Unknown character    (?)	  00000010	         2
---------------------------------------------------


This has been adapated (last two rows have changed) to
allow matches between alignment gaps / ? and any nucleotide:

Nucleotide	    IUPAC code	 Bit-level code	  Value
---------------------------------------------------
A	                  A	      10001000	       136
G	                  G	      01001000	        72
C	                  C	      00101000	        40
T	                  T	      00011000	        24
A or G	              R	      11000000	       192
A or C	              M	      10100000	       160
A or T	              W	      10010000	       144
G or C	              S	      01100000	        96
G or T	              K	      01010000	        80
C or T	              Y	      00110000	        48
A or G or C	          V	      11100000	       224
A or C or T	          H	      10110000	       176
A or G or T	          D	      11010000	       208
G or C or T	          B	      01110000	       112
A or G or C or T      N	      11110000	       240
Alignment gap        (–)	  11110100	       244  ###### CHANGED
Unknown character    (?)	  11110010	       242  ###### CHANGED
---------------------------------------------------


Function	               C code	                      Value returned
------------------------------------------------------------------------------------------
KnownBase(a)	           a & 8	                      8 if a is known surely
IsAdenine(a)	           a == 136	                      1 if a is adenine
IsGuanine(a)	           a == 72	                      1 if a is guanine
IsCytosine(a)       	   a == 40	                      1 if a is cytosine
IsThymine(a)	           a == 24	                      1 if a is thymine
IsPurine(a)	               a & 55	                      0 if a is a purine
IsPyrimidine(a)	           a & 199	                      0 if a is a pyrimidine
DifferentBase(a, b) 	   (a & b) < 16	                  1 if a and b are different surely
SameBase(a, b)	           KnownBase(a) && a == b	      1 if a and b are the same surely
-------------------------------------------------------------------------------------------
=#



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

function make_byte_array()
    byte_array = [0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00, #10
                  0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00, #20
                  0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00, #30
                  0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00, #40
                  0x00,0x00,0x00,0x00,0xf4,0x00,0x00,0x00,0x00,0x00, #50
                  0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00, #60
                  0x00,0x00,0xf2,0x00,0x88,0x70,0x28,0xd0,0x00,0x00, #70
                  0x48,0xb0,0x00,0x00,0x50,0x00,0xa0,0xf0,0x00,0x00, #80
                  0x00,0xc0,0x60,0x18,0x00,0xe0,0x90,0x00,0x30,0x00] #90
    return byte_array
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
