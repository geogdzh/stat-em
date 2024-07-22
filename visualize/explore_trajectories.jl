using CairoMakie, ProgressBars, HDF5, GeoMakie
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl") 

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.

scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
time_history = [x for x in 1850:2014]
time_future = [x for x in 2015:2100]
L1, L2 = 1980, 1032 
l1, l2 = 165, 86
M, N = 192, 96

using_two = true 
non_dim = false  
use_metrics = false
if using_two
    parent_folder = "temp_precip"
else
    parent_folder = "temp"
end
if non_dim
    parent_folder = "nondim"
end
if use_metrics && using_two
    parent_folder = "metrics"
elseif use_metrics && !using_two
    parent_folder = "temp_metrics"
end

#load emulator
d = 10
hfile = h5open("data/$(parent_folder)/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
mean_coefs = read(hfile, "mean_coefs_1")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
close(hfile)

#get gmt list
hfile = h5open("data/ssp245_gmts_50ens.hdf5")
ens_gmt = read(hfile, "ens_gmt")
close(hfile)

gmts = mean(ens_gmt, dims=1)[:]

num_ens_members = 10
trajectories = zeros(M, N, L2, num_ens_members)
for n in 1:num_ens_members
    traj = emulate(gmts, mean_coefs, chol_coefs)
    datatraj = back_to_data(traj, basis) # it's not implemented to work with quadratic yet
    trajectories[:,:,:,n] = shape_data(datatraj[M*N+1:end, :], M, N, true)
end

endtime = time_future[end-59:end]

#load locations
hfile = h5open("data/ground_truth/location_samples_tas_ssp585_50ens.hdf5", "r")
sampling_indices = read(hfile, "indices")
sampling_labels = read(hfile, "labels")
close(hfile)


begin
    fig = Figure()
    for i in 1:3
        for j in 1:2
            ind = (i-1)*2 + j
            ax = Axis(fig[i,j], title=sampling_labels[ind])
            for n in 1:num_ens_members
                shaped = trajectories[:,:,end-59:end,n]
                subset = shaped[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], :]
                avg = mean(subset, dims=(1,2))[:]
                lines!(ax, endtime, avg, color=:black, alpha=0.6)
            end
        end
    end
    # save("data/$(parent_folder)/figs/$(variable)_samples.png", fig)
    display(fig)
end
