using Revise
using FDTDSolver
using LoopVectorization
using BenchmarkTools
using GFlops
using LinuxPerf
using LIKWID

ref_grid = StructuredGrid(LinRange(-20, 15, 10), LinRange(-5, 17, 10), LinRange(-5, 1, 10))

nx, ny, nz = size(ref_grid)


e =  rand(nx, ny, nz, 3)
d = rand(nx, ny, nz, 3)
ϵ =  rand(nx, ny, nz, 3)
b =  rand(nx, ny, nz, 3)
h = rand(nx, ny, nz, 3)
μ =  rand(nx, ny, nz, 3)

j = rand(nx, ny, nz, 3)
σₑ =  rand(nx, ny, nz, 3)

m =  rand(nx, ny, nz, 3)
σₘ =  rand(nx, ny, nz, 3)

cf =  rand(nx, ny, nz, 3)
ca =  rand(nx, ny, nz, 3)

f = rand(nx, ny, nz, 6)


function f()
FaradayStep!(b, e, m, cf, nx, ny, nz)
 #_update_dual_constitutive!(h, b, μ, nx, ny, nz)
 #_update_dual_constitutive!(m, h, σₘ, nx, ny, nz)
 # AmpereStep!(d, h, j, ca, nx, ny, nz)
##_update_primal_constitutive!(j, e, σₑ, nx, ny, nz)
end

@btime f()
#=

cpu = 0 # starts with zero!
LIKWID.PerfMon.init(cpu)
groupid = LIKWID.PerfMon.add_event_set("FLOPS_DP")
LIKWID.PerfMon.setup_counters(groupid)

LIKWID.PerfMon.start_counters()

FaradayStep!(b, e, m,cf, nx, ny, nz)

LIKWID.PerfMon.stop_counters()

mdict = LIKWID.PerfMon.get_metric_results(groupid, cpu)
display(mdict)
println(); flush(stdout);
edict = LIKWID.PerfMon.get_event_results(groupid, cpu)
display(edict)
LIKWID.PerfMon.finalize()
=#