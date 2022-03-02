function _update_dual_constitutive!(h, b, C, nx, ny, nz)
    for k = 2:nz-1
        for j = 2:ny-1
            for i = 2:nx-1
                by = 0.25 * (b.y[i, j, k] + b.y[i-1, j, k] + b.y[i-1, j+1, k] + b.y[i, j+1, k])
                bz = 0.25 * (b.z[i, j, k] + b.z[i, j, k+1] + b.z[i-1, j, k+1] + b.z[i-1, j, k])
                h.x[i, j, k] = C.xx[i, j, k] * b.x[i, j, k] + C.xy[i, j, k] * by + C.xz[i, j, k] * bz


                bx = 0.25 * (b.x[i, j, k] + b.x[i, j-1, k] + b.x[i+1, j-1, k] + b.x[i+1, j, k])
                bz = 0.25 * (b.z[i, j, k] + b.z[i, j-1, k] + b.z[i, j-1, k+1] + b.z[i, j, k+1])

                h.y[i, j, k] = C.yx[i, j, k] * bx + C.yy[i, j, k] * b.y[i, j, k] + C.yz[i, j, k] * bz

                bx = 0.25 * (b.x[i, j, k] + b.x[i+1, j, k] + b.x[i+1, j, k-1] + b.x[i, j, k-1])
                by = 0.25 * (b.y[i, j, k] + b.y[i, j+1, k] + b.y[i, j+1, k-1] + b.y[i, j, k-1])

                h.z[i, j, k] = C.zx[i, j, k] * bx + C.zy[i, j, k] * by + C.zz[i, j, k] * b.z[i, j, k]

            end
        end
    end

    for k = 1:1
        for j = 1:1
            for i = 2:nx-1
                by = 0.25 * (b.y[i, j, k] + b.y[i-1, j, k] + b.y[i-1, j+1, k] + b.y[i, j+1, k])
                bz = 0.25 * (b.z[i, j, k] + b.z[i, j, k+1] + b.z[i-1, j, k+1] + b.z[i-1, j, k])
                h.x[i, j, k] = C.xx[i, j, k] * b.x[i, j, k] + C.xy[i, j, k] * by + C.xz[i, j, k] * bz
            end
        end
    end
    for k = 1:1
        for j = 2:ny-1
            for i = 1:1
                bx = 0.25 * (b.x[i, j, k] + b.x[i, j-1, k] + b.x[i+1, j-1, k] + b.x[i+1, j, k])
                bz = 0.25 * (b.z[i, j, k] + b.z[i, j-1, k] + b.z[i, j-1, k+1] + b.z[i, j, k+1])
                h.y[i, j, k] = C.yx[i, j, k] * bx + C.yy[i, j, k] * b.y[i, j, k] + C.yz[i, j, k] * bz
            end
        end
    end
    for k = 2:ny-1
        for j = 1:1
            for i = 1:1
                bx = 0.25 * (b.x[i, j, k] + b.x[i+1, j, k] + b.x[i+1, j, k-1] + b.x[i, j, k-1])
                by = 0.25 * (b.y[i, j, k] + b.y[i, j+1, k] + b.y[i, j+1, k-1] + b.y[i, j, k-1])

                h.z[i, j, k] = C.zx[i, j, k] * bx + C.zy[i, j, k] * by + C.zz[i, j, k] * b.z[i, j, k]
            end
        end
    end

end