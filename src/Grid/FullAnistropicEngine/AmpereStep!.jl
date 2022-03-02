function AmpereStep!(d, h, J, cb; ref_grid)

    (nx, ny, nz) = size(ref_grid)

    for k = 2:nz-1
        for j = 2:ny-1
            for i = 2:nx-1
                d.x[i, j, k] = d.x[i, j, k] +
                               cb.x[i, j, k] * (h.z[i, j, k] - h.z[i, j-1, k] - h.y[i, j, k] + h.y[i, j, k-1]) - J.x[i, j, k]
                d.y[i, j, k] = d.y[i, j, k] +
                               cb.y[i, j, k] * (h.x[i, j, k] - h.x[i, j, k-1] - h.z[i, j, k] + h.z[i-1, j, k]) - J.y[i, j, k]
                d.z[i, j, k] = d.z[i, j, k] +
                               cb.z[i, j, k] * (h.y[i, j, k] - h.y[i-1, j, k] - h.x[i, j, k] + h.z[i, j-1, k]) - J.z[i, j, k]
            end
        end
    end

    for k = 2:nz-1
        for j = 2:ny-1
            for i = 1:1
                d.x[i, j, k] = d.x[i, j, k] +
                               cb.x[i, j, k] * (h.z[i, j, k] - h.z[i, j-1, k] - h.y[i, j, k] + h.y[i, j, k-1]) - J.x[i, j, k]
            end
        end
    end
    for k = 2:nz-1
        for j = 1:1
            for i = 2:nx-1
                d.y[i, j, k] = d.y[i, j, k] +
                               cb.x[i, j, k] * (h.x[i, j, k] - h.x[i, j, k-1] - h.z[i, j, k] + h.z[i-1, j, k]) - J.y[i, j, k]
            end
        end
    end
    for k = 1:1
        for j = 2:ny-1
            for i = 2:nx-1
                d.z[i, j, k] = d.z[i, j, k] +
                               cb.x[i, j, k] * (h.y[i, j, k] - h.y[i-1, j, k] - h.x[i, j, k] + h.z[i, j-1, k]) - J.z[i, j, k]
            end
        end
    end
end