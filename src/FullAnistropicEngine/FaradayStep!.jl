function FaradayStep!(b, e, m, cb; ref_grid)

    (nx, ny, nz) = size(ref_grid)

    for k = 1:nz-1
        for j = 1:ny-1
            for i = 1:nx-1
                b.x[i, j, k] =
                    b.x[i, j, k] + cb.x[i, j, k] * (e.z[i, j, k] - e.z[i, j+1, k] + e.y[i, j, k+1] - e.y[i, j, k]) - m.x[i, j, k]
                b.y[i, j, k] =
                    b.y[i, j, k] + cb.y[i, j, k] * (e.x[i, j, k] - e.x[i, j, k+1] + e.z[i+1, j, k] - e.z[i, j, k]) - m.y[i, j, k]
                b.z[i, j, k] =
                    b.z[i, j, k] + cb.z[i, j, k] * (e.y[i, j, k] - e.y[i+1, j, k] + e.x[i, j+1, k] - e.x[i, j, k]) -
                    m.z[i, j, k]
            end
        end
    end

    for k = 1:nz-1
        for j = 1:ny-1
            for i = nx:nx
                b.x[i, j, k] =
                    b.x[i, j, k] + cb.x[i, j, k] * (e.z[i, j, k] - e.z[i, j+1, k] + e.y[i, j, k+1] - e.y[i, j, k]) - m.x[i, j, k]
            end
        end
    end
    for k = 1:nz-1
        for j = ny:ny
            for i = 1:nx-1
                b.y[i, j, k] =
                    b.y[i, j, k] + cb.y[i, j, k] * (e.x[i, j, k] - e.x[i, j, k+1] + e.z[i+1, j, k] - e.z[i, j, k]) - m.y[i, j, k]
            end
        end
    end
    for k = nz:nz
        for j = 1:ny-1
            for i = 1:nx-1
                b.z[i, j, k] =
                    b.z[i, j, k] + cb.z[i, j, k] * (e.y[i, j, k] - e.y[i+1, j, k] + e.x[i, j+1, k] - e.x[i, j, k]) -
                    m.z[i, j, k]
            end
        end
    end

end