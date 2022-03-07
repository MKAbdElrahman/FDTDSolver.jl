function AmpereStep!(d, h, J, cb, nx, ny, nz)
    @turbo for k = 2:nz-1
        for j = 2:ny-1
            for i = 2:nx-1
                d[i, j, k, 1] = d[i, j, k, 1] +
                                cb[i, j, k, 1] * (h[i, j, k, 3] - h[i, j-1, k, 3] - h[i, j, k, 2] + h[i, j, k-1, 2]) - J[i, j, k, 1]
                d[i, j, k, 2] = d[i, j, k, 2] +
                                cb[i, j, k, 2] * (h[i, j, k, 1] - h[i, j, k-1, 1] - h[i, j, k, 3] + h[i-1, j, k, 3]) - J[i, j, k, 2]
                d[i, j, k, 3] = d[i, j, k, 3] +
                                cb[i, j, k, 3] * (h[i, j, k, 2] - h[i-1, j, k, 2] - h[i, j, k, 1] + h[i, j-1, k, 1]) - J[i, j, k, 3]
            end
        end
    end
    @turbo for k = 2:nz-1
        for j = 2:ny-1
            for i = 1:1
                d[i, j, k, 1] = d[i, j, k, 1] +
                                cb[i, j, k, 1] * (h[i, j, k, 3] - h[i, j-1, k, 3] - h[i, j, k, 2] + h[i, j, k-1, 2]) - J[i, j, k, 1]
            end
        end
    end
    @turbo for k = 2:nz-1
        for j = 1:1
            for i = 2:nx-1
                d[i, j, k, 2] = d[i, j, k, 2] + cb[i, j, k, 2] * (h[i, j, k, 1] - h[i, j, k-1, 1] - h[i, j, k, 3] + h[i-1, j, k, 3]) - J[i, j, k, 2]
            end
        end
    end
    @turbo for k = 1:1
        for j = 2:ny-1
            for i = 2:nx-1
                d[i, j, k, 3] = d[i, j, k, 3] +
                                cb[i, j, k, 3] * (h[i, j, k, 2] - h[i-1, j, k, 2] - h[i, j, k, 1] + h[i, j-1, k, 1]) - J[i, j, k, 3]
            end
        end
    end
end