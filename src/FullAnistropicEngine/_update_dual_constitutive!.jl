
function _update_dual_constitutive!(H, B, C::Array{T,5}, nx, ny, nz) where T
    @turbo for k = 2:nz-1
        for j = 2:ny-1
            for i = 2:nx-1
                by = 0.25 * (B[i, j, k, 2] + B[i-1, j, k, 2] + B[i-1, j+1, k, 2] + B[i, j+1, k, 2])
                bz = 0.25 * (B[i, j, k, 3] + B[i, j, k+1, 3] + B[i-1, j, k+1, 3] + B[i-1, j, k, 3])
                H[i, j, k, 1] = C[i, j, k, 1, 1] * B[i, j, k, 1] + C[i, j, k, 1, 2] * by + C[i, j, k, 1, 3] * bz

                bx = 0.25 * (B[i, j, k, 1] + B[i, j-1, k, 1] + B[i+1, j-1, k, 1] + B[i+1, j, k, 1])
                bz = 0.25 * (B[i, j, k, 3] + B[i, j-1, k, 3] + B[i, j-1, k+1, 3] + B[i, j, k+1, 3])

                H[i, j, k, 2] = C[i, j, k, 2, 1] * bx + C[i, j, k, 2, 2] * B[i, j, k, 2] + C[i, j, k, 2, 3] * bz

                bx = 0.25 * (B[i, j, k, 1] + B[i+1, j, k, 1] + B[i+1, j, k-1, 1] + B[i, j, k-1, 1])
                by = 0.25 * (B[i, j, k, 2] + B[i, j+1, k, 2] + B[i, j+1, k-1, 2] + B[i, j, k-1, 2])

                H[i, j, k, 3] = C[i, j, k, 3, 1] * bx + C[i, j, k, 3, 2] * by + C[i, j, k, 3, 3] * B[i, j, k, 3]

            end
        end
    end

    @turbo for k = 1:1
        for j = 1:1
            for i = 2:nx-1
                by = 0.25 * (B[i, j, k, 2] + B[i-1, j, k, 2] + B[i-1, j+1, k, 2] + B[i, j+1, k, 2])
                bz = 0.25 * (B[i, j, k, 3] + B[i, j, k+1, 3] + B[i-1, j, k+1, 3] + B[i-1, j, k, 3])
                H[i, j, k, 1] = C[i, j, k, 1, 1] * B[i, j, k, 1] + C[i, j, k, 1, 2] * by + C[i, j, k, 1, 3] * bz

            end
        end
    end
    @turbo for k = 1:1
        for j = 2:ny-1
            for i = 1:1
                bx = 0.25 * (B[i, j, k, 1] + B[i, j-1, k, 1] + B[i+1, j-1, k, 1] + B[i+1, j, k, 1])
                bz = 0.25 * (B[i, j, k, 3] + B[i, j-1, k, 3] + B[i, j-1, k+1, 3] + B[i, j, k+1, 3])
                H[i, j, k, 2] = C[i, j, k, 2, 1] * bx + C[i, j, k, 2, 2] * B[i, j, k, 2] + C[i, j, k, 2, 3] * bz
            end
        end
    end
    @turbo for k = 2:nz-1
        for j = 1:1
            for i = 1:1
                bx = 0.25 * (B[i, j, k, 1] + B[i+1, j, k, 1] + B[i+1, j, k-1, 1] + B[i, j, k-1, 1])
                by = 0.25 * (B[i, j, k, 2] + B[i, j+1, k, 2] + B[i, j+1, k-1, 2] + B[i, j, k-1, 2])
                H[i, j, k, 3] = C[i, j, k, 3, 1] * bx + C[i, j, k, 3, 2] * by + C[i, j, k, 3, 3] * B[i, j, k, 3]
            end
        end
    end

end



function _update_dual_constitutive!(H, B, C::Array{T,4}, nx, ny, nz) where T
    @turbo for p in 1:3
        for k = 2:nz-1
            for j = 2:ny-1
                for i = 2:nx-1
                    H[i, j, k, p] = C[i, j, k, p] * B[i, j, k, p]
                end
            end
        end
    end

    @turbo for k = 1:1
        for j = 1:1
            for i = 2:nx-1
                H[i, j, k, 1] = C[i, j, k, 1] * B[i, j, k, 1]
            end
        end
    end
    @turbo for k = 1:1
        for j = 2:ny-1
            for i = 1:1
                H[i, j, k, 2] = C[i, j, k, 2] * B[i, j, k, 2]
            end
        end
    end
    @turbo for k = 2:nz-1
        for j = 1:1
            for i = 1:1
                H[i, j, k, 3] = C[i, j, k, 3] * B[i, j, k, 3]
            end
        end
    end
end