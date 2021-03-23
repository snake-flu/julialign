include("functions_collapse.jl")

function pairsnp_square(nuc_bit_array_A, nuc_bit_V)

    diff_list_A = get_list_of_differences(nuc_bit_array_A, nuc_bit_V)

    snp_count_array = zeros(Int, size(nuc_bit_array_A, 2), size(nuc_bit_array_A, 2))

    # For parallelising:
    n_Chunks = Threads.nthreads()
    if n_Chunks > size(nuc_bit_array_A, 2) / 2
        n_Chunks::Int = floor(size(nuc_bit_array_A, 2) / 2)
    end
    ranges = get_ranges(size(nuc_bit_array_A, 2), n_Chunks)

    Threads.@threads for range in ranges
        for a in range[1]:range[2]
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
    end

    return snp_count_array
end

function pairsnp_rectangle(nuc_bit_array_A, nuc_bit_array_B, nuc_bit_V)

    diff_list_A = get_list_of_differences(nuc_bit_array_A, nuc_bit_V)
    diff_list_B = get_list_of_differences(nuc_bit_array_B, nuc_bit_V)

    snp_count_array = zeros(Int, size(nuc_bit_array_B, 2), size(nuc_bit_array_A, 2))

    Threads.@threads for a in 1:size(nuc_bit_array_A, 2)
        for b in 1:size(nuc_bit_array_B, 2)

            snp_counter = 0

            for r in diff_list_A[a]
                @inbounds x = nuc_bit_array_A[r,a]
                @inbounds y = nuc_bit_array_B[r,b]
                different = (x & y) < 16

                if different
                    snp_counter+=1
                end
            end

            for r in diff_list_B[b]
                if r in diff_list_A[a]
                    continue
                end

                @inbounds x = nuc_bit_array_A[r,a]
                @inbounds y = nuc_bit_array_B[r,b]
                different = (x & y) < 16

                if different
                    snp_counter+=1
                end
            end

            snp_count_array[b,a] = snp_counter
        end
    end

    return snp_count_array
end

function write_pairsnp_square(outfile, int_array, seq_name_vector)

    open(outfile, "w") do io
        @inbounds for i in 1:size(int_array, 2)
            write(io, seq_name_vector[i])
            write(io, ",")

            write(io, join(int_array[:,i], ",") * '\n')
        end
    end
end

function write_pairsnp_rectangle(outfile, int_array, seq_name_vector_A, seq_name_vector_B)

    open(outfile, "w") do io
        write(io, ",")
        write(io, join(seq_name_vector_B, ",") * "\n")

        @inbounds for i in 1:size(int_array, 2)
            write(io, seq_name_vector_A[i])
            write(io, ",")
            write(io, join(int_array[:,i], ",") * '\n')
        end
    end
end
