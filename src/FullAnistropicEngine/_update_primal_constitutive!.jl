function _update_primal_constitutive!(e, d, C::Array{T,5}, nx, ny, nz) where {T}
    @turbo for k = 2:nz-1
        for j = 2:ny-1
            for i = 2:nx-1
                dy = 0.25 * (d[i, j, k, 2] + d[i+1, j, k, 2] + d[i+1, j-1, k, 2] + d[i, j-1, k, 2])
                dz = 0.25 * (d[i, j, k, 3] + d[i+1, j, k, 3] + d[i+1, j, k-1, 3] + d[i, j, k-1, 3])
                e[i, j, k, 1] = C[i, j, k, 1, 1] * d[i, j, k, 1] + C[i, j, k, 1, 2] * dy + C[i, j, k, 1, 3] * dz


                dx = 0.25 * (d[i, j, k, 1] + d[i, j+1, k, 1] + d[i-1, j, k, 1] + d[i-1, j+1, k, 1])
                dz = 0.25 * (d[i, j, k, 3] + d[i+1, j, k, 3] + d[i+1, j, k-1, 3] + d[i, j, k-1, 3])

                e[i, j, k, 2] = C[i, j, k, 2, 1] * dx + C[i, j, k, 2, 2] * d[i, j, k, 2] + C[i, j, k, 2, 3] * dz

                dx = 0.25 * (d[i, j, k, 1] + d[i, j, k+1, 1] + d[i-1, j, k+1, 1] + d[i-1, j, k, 1])
                dy = 0.25 * (d[i, j, k, 2] + d[i, j, k+1, 2] + d[i, j-1, k+1, 2] + d[i, j-1, k, 2])

                e[i, j, k, 3] = C[i, j, k, 3, 1] * dx + C[i, j, k, 3, 2] * dy + C[i, j, k, 3, 3] * d[i, j, k, 3]

            end
        end
    end

    @turbo for k = 2:nz-1
        for j = 2:ny-1
            for i = 1:1

                dy = 0.25 * (d[i, j, k, 2] + d[i+1, j, k, 2] + d[i+1, j-1, k, 2] + d[i, j-1, k, 2])
                dz = 0.25 * (d[i, j, k, 3] + d[i+1, j, k, 3] + d[i+1, j, k-1, 3] + d[i, j, k-1, 3])
                e[i, j, k, 1] = C[i, j, k, 1, 1] * d[i, j, k, 1] + C[i, j, k, 1, 2] * dy + C[i, j, k, 1, 3] * dz

            end
        end
    end
    @turbo for k = 2:nz-1
        for j = 1:1
            for i = 2:nx-1

                dx = 0.25 * (d[i, j, k, 1] + d[i, j+1, k, 1] + d[i-1, j, k, 1] + d[i-1, j+1, k, 1])
                dz = 0.25 * (d[i, j, k, 3] + d[i+1, j, k, 3] + d[i+1, j, k-1, 3] + d[i, j, k-1, 3])
                e[i, j, k, 2] = C[i, j, k, 2, 1] * dx + C[i, j, k, 2, 2] * d[i, j, k, 2] + C[i, j, k, 2, 3] * dz

            end
        end
    end
    @turbo for k = 1:1
        for j = 2:ny-1
            for i = 2:nx-1
                dx = 0.25 * (d[i, j, k, 1] + d[i, j, k+1, 1] + d[i-1, j, k+1, 1] + d[i-1, j, k, 1])
                dy = 0.25 * (d[i, j, k, 2] + d[i, j, k+1, 2] + d[i, j-1, k+1, 2] + d[i, j-1, k, 2])
                e[i, j, k, 3] = C[i, j, k, 3, 1] * dx + C[i, j, k, 3, 2] * dy + C[i, j, k, 3, 3] * d[i, j, k, 3]

            end
        end
    end
end

function _update_primal_constitutive!(e, d, C::Array{T,4}, nx, ny, nz) where {T}
    @turbo for p in 1:3
        for k = 2:nz-1
            for j = 2:ny-1
                for i = 2:nx-1
                    e[i, j, k, p] = C[i, j, k, p] * d[i, j, k, p]
                end
            end
        end
    end
    @turbo for k = 2:nz-1
        for j = 2:ny-1
            for i = 1:1
                e[i, j, k, 1] = C[i, j, k, 1] * d[i, j, k, 1]
            end
        end
    end
    @turbo for k = 2:nz-1
        for j = 1:1
            for i = 2:nx-1
                e[i, j, k, 2] = C[i, j, k, 2] * d[i, j, k, 2]
            end
        end
    end
    @turbo for k = 1:1
        for j = 2:ny-1
            for i = 2:nx-1
                e[i, j, k, 3] = C[i, j, k, 3] * d[i, j, k, 3]
            end
        end
    end
end