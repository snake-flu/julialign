include("functions_align.jl")
include("functions_encoding.jl")

function is_column_biallelic(nuc_bit_vector)

    nucs = Set{UInt8}()
    counter = 0

    for nuc in nuc_bit_vector
        if nuc & 8 != 8
            continue
        end
        if ! (nuc in nucs)
            push!(nucs, nuc)
            counter += 1
            if counter > 2
                return false
            end
        end
    end

    if counter == 2
        return true
    end

    return false
end

function collapse_identical_haplotypes(nuc_bit_array)

    hap_dict = Dict{Vector{UInt8}, Vector{Int}}()

    for i in 1:size(nuc_bit_array, 1)
        if in(nuc_bit_array[i,:], keys(hap_dict))
            push!(hap_dict[nuc_bit_array[i,:]], i)
        else
            hap_dict[nuc_bit_array[i,:]] = [i]
        end
    end

    return hap_dict
end

function four_haplotypes_questionmark(hap_array)
    s = Set()
    for i in 1:size(hap_array, 1)
        if ! (hap_array[i,1] & 8 == 8)
            continue
        end
        if ! (hap_array[i,2] & 8 == 8)
            continue
        end
        push!(s, hap_array[i,:])
    end
    if length(s) == 4
        return true
    else
        return false
    end
end

# reduce / get an array of biallelic ATGC sites only
function get_biallelic_columns(nuc_bit_array, fasta_IDs)

    nuc_bit_array = transpose(nuc_bit_array)

    hap_dict = collapse_identical_haplotypes(nuc_bit_array)

    rows = Vector{Int}()
    for v in values(hap_dict)
        push!(rows, v[1])
    end

    if length(rows) < 4
        e = error("not enough unique haplotypes to perform the 4-gamete test")
        throw(e)
    end

    sort!(rows)

    unique_nuc_bit_array = nuc_bit_array[rows, :]

    bi_idx = Array{Int, 1}(undef, 0)

    for i in 1:size(unique_nuc_bit_array, 2)
        if is_column_biallelic(unique_nuc_bit_array[:,i])
            push!(bi_idx, i)
        end
    end

    println(length(bi_idx))
    println()

    bi_view = view(unique_nuc_bit_array, :, bi_idx)

    nd = make_nuc_dict()

    for i in 1:size(bi_view, 1)
        for j in 1:size(bi_view, 2)
            print(nd[bi_view[i,j]])
        end
        print(" " * fasta_IDs[rows[i]])
        println()
    end

    println()

    D = zeros(Int8, length(bi_idx), length(bi_idx))

    state = 0
    for i in 1:(size(bi_view, 2) - 1)
        for j in (i + 1):size(bi_view, 2)
            gametes = view(bi_view, :, [i,j])
            four_gamete_test = four_haplotypes_questionmark(gametes)
            if four_gamete_test
                D[i,j] = 1
                D[j,i] = 1
            end
        end
    end

    for i in 1:size(D, 1)
        for j in 1:size(D, 2)
            print(D[i,j])
        end
        print(" ")
        print(bi_idx[i])
        println()
    end

end





















#
