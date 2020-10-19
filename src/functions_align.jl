include("functions_encoding.jl")


# get_index_of_diff_sites returns an array of integers which are the sites
# at which two vectors of bit-level encoded nucleotides differ
# NB nuc_bit_V_A should be the reference, as it stands, because we ignore N/-/?
# in nuc_bit_V_B
function get_index_of_diff_sites(nuc_bit_V_A, nuc_bit_V_B)

    A = Array{Int, 1}()

    for i in 1:length(nuc_bit_V_A)

        @inbounds x = nuc_bit_V_A[i]
        @inbounds y = nuc_bit_V_B[i]

        if x != y && y < 240
            push!(A, i)
        end
    end

    return A
end


# get_list_of_differences applies get_index_of_diff_sites between nuc_bit_V_A
# and every sequence stored in nuc_bit_array
function get_list_of_differences(nuc_bit_array, nuc_bit_V_A)

    L = Array{Array{Int, 1}, 1}(undef, size(nuc_bit_array, 2))

    for i in 1:length(L)
        L[i] = get_index_of_diff_sites(nuc_bit_V_A, view(nuc_bit_array, :, i))
    end

    return L
end

struct fasta_record
    id::String
    description::String
    seq::String
    n::Int
end

function get_alignment_dimensions(filepath::AbstractString)

    n = 0 # number of sequences
    l = 0 # length of sequence
    first = true

    open(filepath, "r") do io
        while !eof(io)

            if n > 1
                first = false
            end

            line = readline(io)

            if string(line[1]) == ">"
                n+=1
            end

            if first && string(line[1]) != ">"
                l+=length(line)
            end
        end

    end

    return l, n
end

function get_fasta_descriptions(filepath::AbstractString)

    A = Array{String,1}()

    open(filepath, "r") do io
        while !eof(io)
            l = readline(io)
            if string(l[1]) == ">"
                description = lstrip(l, '>')
                id = split(description)[1]
                push!(A, description)
            end
        end
    close(io)
    end

    return A
end

function read_fasta_alignment(filepath::AbstractString, channel::Channel)

    open(filepath, "r") do io

        n = 0

        seq_buffer = String("")
        id = String("")
        description = String("")

        while !eof(io)

            l = readline(io)
            n += 1

            if n == 1
                if string(l[1]) != ">"
                    throw("badly formed fasta file")
                end

                description = lstrip(l, '>')
                id = split(description)[1]

                first = false

            elseif eof(io)

                seq_buffer = seq_buffer * uppercase(l)

                fr = fasta_record(id, description, seq_buffer, n)
                put!(channel, fr)

                break

            elseif string(l[1]) == ">"

                fr = fasta_record(id, description, seq_buffer, n)
                put!(channel, fr)

                description = lstrip(l, '>')
                id = split(description)[1]
                seq_buffer = String("")

            else

                seq_buffer = seq_buffer * uppercase(l)

            end
        end

    close(io)
    end
end

function populate_byte_array(filepath::AbstractString)
    byte_dict = make_byte_dict()

    alignment_dim = get_alignment_dimensions(filepath)
    height = alignment_dim[1]
    width = alignment_dim[2]

    A = Array{UInt8,2}(undef, height, width)

    ch = Channel{fasta_record}((channel_arg) -> read_fasta_alignment(filepath, channel_arg))


    j = 1
    for record in ch
        sequence = record.seq
        i = 1
        for nuc in sequence
            # byte = byte_dict[nuc]
            A[i, j] = byte_dict[nuc]
            i+=1
        end
        j+=1
    end

    return A
end

function get_seq_from_1D_byte_array(byte_V)
    nuc_dict = make_nuc_dict()

    char_V = Array{Char, 1}(undef, length(byte_V))

    for i in 1:size(byte_V, 1)
        char_V[i] = nuc_dict[byte_V[i]]
    end

    return join(char_V)
end

function score_alignment_column(nuc_bit_sequence)
    score_dict = make_score_dict()

    score = 0::Int64
    for nuc in nuc_bit_sequence
        score+=score_dict[nuc]
    end
    return score
end

function score_alignment(nuc_bit_array)

    A = fill(0::Int64, size(nuc_bit_array, 2))

    for i in 1:size(nuc_bit_array, 2)
        A[i] = score_alignment_column(view(nuc_bit_array, :, i))
    end

    return A
end
