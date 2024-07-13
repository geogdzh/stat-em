using CairoMakie, ProgressBars, HDF5, GeoMakie
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl") 

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.

filepr = file_head*"ssp585/pr/r1i1p1f1_ssp585_pr.nc"
pr = ncData(filepr, "pr")
# prslice = pr(208, 210, 62, 65)
prslice = pr(285, 289 , 40, 44)
prdata = prslice.data .* 86400 #convert to mm/day

size(prdata)
avg = mean(log.(prdata), dims=(1,2))[:]
hist(avg, bins=20, xlabel="Precipitation (mm/day)", ylabel="Frequency")