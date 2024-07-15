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


### let's identify the locations
ts3(7,9,11,12).lonvec

function get_indices(lon1, lon2, lat1, lat2, ts)
    lons = ts(lon1, lon2, lat1, lat2).lonvec
    lats = ts(lon1, lon2, lat1, lat2).latvec
    return [findfirst(x -> x == lons[1], ts.lonvec), 
    findfirst(x -> x == lons[end], ts.lonvec), 
    findfirst(x -> x == lats[1], ts.latvec), 
    findfirst(x -> x == lats[end], ts.latvec)]
end

sampling_indices = zeros(Int64, (6, 4))
sampling_indices[1,:] = get_indices(272,273,41,42, ts3) #chicago
sampling_indices[2,:] = get_indices(300,301,-4,-3, ts3) #manaus, brazil
sampling_indices[3,:] = get_indices(7,9,11,12, ts3) #nigeria
sampling_indices[4,:] = get_indices(26,28,68,70, ts3) #finland
sampling_indices[5,:] = get_indices(90,92,23,25, ts3) #dhaka, bhangladesh
sampling_indices[6,:] = get_indices(174, 176, -38, -36, ts3) #auckland, nz

sampling_labels = ["Chicago, USA", "Manaus, Brazil", "Kano, Nigeria", "Inari, Finland", "Dhaka, Bhangladesh", "Auckland, NZ"]

#### start with: real data PDFs for precip and temp
# num_ens_members = 50
# for variable in ["tas", "pr"]
#     points = zeros(6, num_ens_members, L2)
#     for n in 1:num_ens_members
#         file = file_head*"ssp585/$(variable)/r$(n)i1p1f1_ssp585_$(variable).nc"
#         ts = ncData(file, variable)
#         data = variable == "pr" ? log.(ts.data .* 86400) : ts.data .* 1.0
#         for ind in 1:6
#             subset = data[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], :]
#             avg = mean(subset, dims=(1,2))[:]
#             points[ind, n, :] = avg
#         end
#     end
#     hfile = h5open("data/ground_truth/location_samples_$(variable)_ssp585_50ens.hdf5", "w")
#     write(hfile, "points", points)
#     write(hfile, "indices", sampling_indices)
#     write(hfile, "labels", sampling_labels)
#     close(hfile)
# end

#####


variable = "tas"

hfile = h5open("data/ground_truth/vars_$(variable)_ssp585_50ens.hdf5", "r") # true CMIP vars incl means
true_var = read(hfile, "true_var")
true_ens_mean = read(hfile, "true_ens_mean")
close(hfile)
sample_means = zeros(6, l2)
sample_vars = zeros(6, l2)
for ind in 1:6
    slice = true_ens_mean[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], 1:12:end] #only jan
    sample_means[ind, :] = mean(slice, dims=(1,2))[:]
    slice = true_var[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], 1:12:end] #only jan
    sample_vars[ind, :] = mean(slice, dims=(1,2))[:]    #is this a valid way of doing it??
end


hfile = h5open("data/ground_truth/location_samples_$(variable)_ssp585_50ens.hdf5", "r")
points = read(hfile, "points")
close(hfile)
begin
    data = points[:,:,1:12:end]# looking only at January!
    fig = Figure()
    for i in 1:3
        for j in 1:2
            ind = (i-1)*2 + j
            ax = Axis(fig[i,j], title=sampling_labels[ind])
            samples = reshape(data[ind,:,:], (num_ens_members*l2,))
            hist!(ax, samples, bins=20, xlabel="Temperature (K)", ylabel="Frequency", normalization=:pdf)

            # add the distribution
            μ = mean(sample_means[ind, :]) #is this how it works???
            var = mean(sample_vars[ind, :])
            dist = Normal(μ, sqrt(var))
            plot!(ax, dist, label="True distribution", color=:red)
        end
    end
    display(fig)
end







####



hfile = h5open("data/ground_truth/vars_tas_ssp585_50ens.hdf5", "r") # true CMIP vars
true_var = read(hfile, "true_var")
true_ens_mean = read(hfile, "true_ens_mean")
close(hfile)

hfile = h5open("../cmip-mc/data/temp_precip/vars_temp_ssp585_50ens.hdf5", "r") # true CMIP vars
true_var_og = read(hfile, "true_var")
true_ens_mean_og = read(hfile, "true_ens_mean")
close(hfile)

true_var == true_var_og
true_ens_mean == true_ens_mean_og

##

hfile = h5open("data/temp_precip/ens_vars/ens_vars_ssp585_10d.hdf5", "r")
ens_vars_tas_10 = read(hfile, "ens_vars_tas_10")
ens_means_tas_10 = read(hfile, "ens_means_tas_10")
close(hfile)

hfile = h5open("../cmip-mc/data/temp_precip/ens_vars_withpr_ssp585.hdf5", "r") # true CMIP vars
ens_vars_tas_20_og = read(hfile, "ens_vars_tas_20")
ens_means_tas_20_og = read(hfile, "ens_means_tas_20")
close(hfile)

diff = ens_vars_tas_10 .- ens_vars_tas_20_og
for i in 1:100
    fig = Figure()
    ax = Axis(fig[1,1])
    heatmap!(ax,diff[:,:,i], colormap=:coolwarm, clims=(-1, 1))
    display(fig)
end