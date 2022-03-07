function FaradayStep!(b, e, m, cb, nx, ny, nz)
    @turbo for k = 1:nz-1
        for j = 1:ny-1
            for i = 1:nx-1
                b[i, j, k, 1] =
                    b[i, j, k, 1] + cb[i, j, k, 1] * (e[i, j, k, 3] - e[i, j+1, k, 3] + e[i, j, k+1, 2] - e[i, j, k, 2]) - m[i, j, k, 1]
                b[i, j, k, 2] =
                    b[i, j, k, 2] + cb[i, j, k, 2] * (e[i, j, k, 1] - e[i, j, k+1, 1] + e[i+1, j, k, 3] - e[i, j, k, 3]) - m[i, j, k, 2]
                b[i, j, k, 3] =
                    b[i, j, k, 3] + cb[i, j, k, 3] * (e[i, j, k, 2] - e[i+1, j, k, 2] + e[i, j+1, k, 1] - e[i, j, k, 1]) - m[i, j, k, 3]
            end
        end
    end
    @turbo for k = 1:nz-1
        for j = 1:ny-1
            for i = nx:nx
                b[i, j, k, 1] =
                    b[i, j, k, 1] + cb[i, j, k, 1] * (e[i, j, k, 3] - e[i, j+1, k, 3] + e[i, j, k+1, 2] - e[i, j, k, 2]) - m[i, j, k, 1]
            end
        end
    end
    @turbo for k = 1:nz-1
        for j = ny:ny
            for i = 1:nx-1
                b[i, j, k, 2] =
                    b[i, j, k, 2] + cb[i, j, k, 2] * (e[i, j, k, 1] - e[i, j, k+1, 1] + e[i+1, j, k, 3] - e[i, j, k, 3]) - m[i, j, k, 2]
            end
        end
    end
    @turbo for k = nz:nz
        for j = 1:ny-1
            for i = 1:nx-1
                b[i, j, k, 3] =
                    b[i, j, k, 3] + cb[i, j, k, 3] * (e[i, j, k, 2] - e[i+1, j, k, 2] + e[i, j+1, k, 1] - e[i, j, k, 1]) - m[i, j, k, 3]
            end
        end
    end
end



            #    ex hx ey hy ez hz 
function FaradayStep!(f, m, cb, nx, ny, nz)
    @turbo for k = 1:nz-1
        for j = 1:ny-1
            for i = 1:nx-1
                f[i, j, k, 2] =
                    f[i, j, k, 2] + cb[i, j, k, 1] * (f[i, j, k, 5] - f[i, j+1, k, 5] + f[i, j, k+1, 3] - f[i, j, k, 3]) - m[i, j, k, 1]
                f[i, j, k, 4] =
                    f[i, j, k, 4] + cb[i, j, k, 2] * (f[i, j, k, 1] - f[i, j, k+1, 1] + f[i+1, j, k, 5] - f[i, j, k, 5]) - m[i, j, k, 2]
                f[i, j, k, 6] =
                    f[i, j, k, 6] + cb[i, j, k, 3] * (f[i, j, k, 3] - f[i+1, j, k, 3] + f[i, j+1, k, 1] - f[i, j, k, 1]) - m[i, j, k, 3]
            end
        end
    end
    @turbo for k = 1:nz-1
        for j = 1:ny-1
            for i = nx:nx
                f[i, j, k, 4] =
                f[i, j, k, 4] + cb[i, j, k, 1] * (f[i, j, k, 3] - f[i, j+1, k, 3] + f[i, j, k+1, 2] - f[i, j, k, 2]) - m[i, j, k, 1]
            end
        end
    end
    @turbo for k = 1:nz-1
        for j = ny:ny
            for i = 1:nx-1
                f[i, j, k, 5] =
                f[i, j, k, 5] + cb[i, j, k, 2] * (f[i, j, k, 1] - f[i, j, k+1, 1] + f[i+1, j, k, 3] - f[i, j, k, 3]) - m[i, j, k, 2]     
            end
        end
    end
    @turbo for k = nz:nz
        for j = 1:ny-1
            for i = 1:nx-1
                f[i, j, k, 6] =
                f[i, j, k, 6] + cb[i, j, k, 3] * (f[i, j, k, 2] - f[i+1, j, k, 2] + f[i, j+1, k, 1] - f[i, j, k, 1]) - m[i, j, k, 3]
            end
        end
    end
end