

# function pairsnp(nuc_bit_array_A, nuc_bit_array_B, nuc_bit_V)
function pairsnp(nuc_bit_array_A, nuc_bit_V)

    diff_list_A = get_list_of_differences(nuc_bit_array_A, nuc_bit_V)

    snp_count_array = zeros(Int, size(nuc_bit_array_A, 2), size(nuc_bit_array_A, 2))

    for a in 1:size(nuc_bit_array_A, 2) - 1

        for b in (a + 1):size(nuc_bit_array_A, 2)

            snp_counter = 0

            for r in diff_list_A[a]
                @inbounds x = nuc_bit_array_A[r,a]
                @inbounds y = nuc_bit_array_A[r,b]
                different = (x & y) < 16

                if different
                    snp_counter+=1
                end
            end

            for r in diff_list_A[b]
                if r in diff_list_A[a]
                    continue
                end

                @inbounds x = nuc_bit_array_A[r,a]
                @inbounds y = nuc_bit_array_A[r,b]
                different = (x & y) < 16

                if different
                    snp_counter+=1
                end
            end

            snp_count_array[a,b] = snp_counter
            snp_count_array[b,a] = snp_counter
        end
    end

    return snp_count_array
end

function write_pairsnp(outfile, int_array, seq_name_vector)

    open(outfile, "w") do io
        @inbounds for i in 1:size(int_array, 2)
            write(io, seq_name_vector[i])
            write(io, ",")

            write(io, join(int_array[:,i], ",") * '\n')
        end
    end
end
