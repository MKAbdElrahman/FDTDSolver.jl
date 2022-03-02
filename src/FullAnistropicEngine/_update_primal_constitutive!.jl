function _update_primal_constitutive!(e, d, C, nx, ny, nz)
    for k = 2:nz-1
        for j = 2:ny-1
            for i = 2:nx-1
                dy = 0.25 * (d.y[i, j, k] + d.y[i+1, j, k] + d.y[i+1, j-1, k] + d.y[i, j-1, k])
                dz = 0.25 * (d.z[i, j, k] + d.z[i+1, j, k] + d.z[i+1, j, k-1] + d.z[i, j, k-1])
                e.x[i, j, k] = C.xx[i, j, k] * d.x[i, j, k] + C.xy[i, j, k] * dy + C.xz[i, j, k] * dz


                dx = 0.25 * (d.x[i, j, k] + d.x[i, j+1, k] + d.x[i-1, j, k] + d.x[i-1, j+1, k])
                dz = 0.25 * (d.z[i, j, k] + d.z[i+1, j, k] + d.z[i+1, j, k-1] + d.z[i, j, k-1])

                e.y[i, j, k] = C.yx[i, j, k] * dx + C.yy[i, j, k] * d.y[i, j, k] + C.yz[i, j, k] * dz

                dx = 0.25 * (d.x[i, j, k] + d.x[i, j, k+1] + d.x[i-1, j, k+1] + d.x[i-1, j, k])
                dy = 0.25 * (d.y[i, j, k] + d.y[i, j, k+1] + d.y[i, j-1, k+1] + d.y[i, j-1, k])

                e.z[i, j, k] = C.zx[i, j, k] * dx + C.zy[i, j, k] * dy + C.yz[i, j, k] * d.z[i, j, k]

            end
        end
    end

    for k = 2:nz-1
        for j = 2:ny-1
            for i = 1:1

                dy = 0.25 * (d.y[i, j, k] + d.y[i+1, j, k] + d.y[i+1, j-1, k] + d.y[i, j-1, k])
                dz = 0.25 * (d.z[i, j, k] + d.z[i+1, j, k] + d.z[i+1, j, k-1] + d.z[i, j, k-1])
                e.x[i, j, k] = C.xx[i, j, k] * d.x[i, j, k] + C.xy[i, j, k] * dy + C.xz[i, j, k] * dz

            end
        end
    end
    for k = 2:nz-1
        for j = 1:1
            for i = 2:nx-1

                dx = 0.25 * (d.x[i, j, k] + d.x[i, j+1, k] + d.x[i-1, j, k] + d.x[i-1, j+1, k])
                dz = 0.25 * (d.z[i, j, k] + d.z[i, j+1, k] + d.z[i, j+1, k-1] + d.z[i, j, k-1])

                e.y[i, j, k] = C.yx[i, j, k] * dx + C.yy[i, j, k] * d.y[i, j, k] + C.yz[i, j, k] * dz
            end
        end
    end
    for k = 1:1
        for j = 2:ny-1
            for i = 2:nx-1
                dx = 0.25 * (d.x[i, j, k] + d.x[i, j, k+1] + d.x[i-1, j, k+1] + d.x[i-1, j, k])
                dy = 0.25 * (d.y[i, j, k] + d.y[i, j, k+1] + d.y[i, j-1, k+1] + d.y[i, j-1, k])

                e.z[i, j, k] = C.zx[i, j, k] * dx + C.zy[i, j, k] * dy + C.yz[i, j, k] * d.z[i, j, k]
            end
        end
    end
end
