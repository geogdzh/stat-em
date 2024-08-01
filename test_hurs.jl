using CairoMakie, ProgressBars, HDF5, GeoMakie
include("utils/data_util.jl")
include("utils/eof_util.jl")
include("utils/emulator_util.jl") 

parent_folder = "hurs_only"

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"

file_tail = "historical/hurs/r10i1p1f1_historical_hurs.nc"
ts = ncData(file_head*file_tail, "hurs")
M, N, L1 = size(ts.data)

X = reshape_data(ts.data)

file_tail = "ssp585/hurs/r10i1p1f1_ssp585_hurs.nc"
ts2 = ncData(file_head*file_tail, "hurs")
X2 = reshape_data(ts2.data)
M, N, L2 = size(ts2.data)

fullX = hcat(X, X2)

U, S, V = svd(fullX)
d = 2000
basis = U[:, 1:d]

hfile = h5open("data/$(parent_folder)/basis_$(d)d.hdf5", "w") #but to use a smaller basis, can just take fewer modes
write(hfile, "basis", basis)
close(hfile)

scenarios = ["historical", "ssp585"]

d = 20
basis = basis[:, 1:d]

ss = sum(S.^2)
cum_var = cumsum(S.^2 ./ ss)

begin
    fig = Figure(resolution=(600, 400))
    ax = Axis(fig[1,1], xlabel="Mode", ylabel="Cumulative variance explained", title="Cumulative variance explained by the modes")
    lines!(ax, cum_var, color=:black, linestyle=:solid)
    display(fig)
    # save("figs/cum_var.png", fig)
end

###

variable = "hurs"

label = variable == "tas" ? variable : "two"

# get the ground truth
hfile = h5open("data/ground_truth/location_samples_$(variable)_ssp585_$(num_ens_members)ens.hdf5", "r")
sampling_indices = read(hfile, "sampling_indices")
sampling_labels = read(hfile, "sampling_labels")
end_points = read(hfile, "end_points")
sample_means = read(hfile, "sample_means")
sample_vars = read(hfile, "sample_vars")
close(hfile)


#### start with: real data PDFs for last ten years of SSP585
end_points = zeros(6, num_ens_members, 120)
for n in 1:num_ens_members
    file = file_head*"ssp585/$(variable)/r$(n)i1p1f1_ssp585_$(variable).nc"
    ts = ncData(file, variable)
    data = apply_transform(ts.data, variable; hurs_option=hurs_option)
    for ind in 1:6
        subset = data[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], end-119:end]
        avg = mean(subset, dims=(1,2))[:]
        end_points[ind, n, :] = avg
    end
end
# wfile = h5open("data/ground_truth/location_samples_$(variable)_ssp585_$(num_ens_members)ens.hdf5", "w")
# write(wfile, "end_points", end_points)
# write(wfile, "sampling_indices", sampling_indices)
# write(wfile, "sampling_labels", sampling_labels)

# normal approximation to ground truth
hfile = h5open("data/ground_truth/vars_$(variable)_ssp585_$(num_ens_members)ens.hdf5", "r") # true CMIP vars incl means
true_var = read(hfile, "true_var")
true_ens_mean = read(hfile, "true_ens_mean")
close(hfile)
sample_means = zeros(6, l2)
sample_vars = zeros(6, l2)
for ind in 1:6
    slice = true_ens_mean[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], 1:12:end] #only jan!!
    sample_means[ind, :] = mean(slice, dims=(1,2))[:]
    slice = true_var[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], 1:12:end] #only jan!!
    sample_vars[ind, :] = mean(slice, dims=(1,2))[:] 
end
# write(wfile, "sample_means", sample_means)
# write(wfile, "sample_vars", sample_vars)
# close(wfile)

########################

# hfile = h5open("data/ssp585_gmts_28ens.hdf5", "r")
# ens_gmt = read(hfile, "ens_gmt")
# close(hfile)
# ens_gmt = mean(ens_gmt, dims=1)[:]

parent_folder = "temp_hurs"
d = 100
scenario = "ssp585"
num_ens_members = 28

hfile = h5open("data/$(parent_folder)/training_data_ssp585_$(d)d_$(num_ens_members)ens.hdf5", "r")
ens_gmt = read(hfile, "ens_gmt")
close(hfile)
ens_gmt = mean(ens_gmt, dims=1)[:]

chol_coefs, chols = get_chol_coefs(ens_gmt, "$("ssp585")_$("100")d"; return_chols=true, offload=false, parent_folder=nothing)

###

hfile = h5open("data/process/chols_$(scenario)_$(d)d.hdf5", "r")
chols = read(hfile, "chols")
chol_coefs = read(hfile, "chol_coefs")
ens_gmt = read(hfile, "ens_gmt")
close(hfile)
ens_gmt = ens_gmt[:]

begin
    shift  = 90
    fig = Figure(resolution=(1000, 1000))
    for i in 1:5
        for j in 1:5
            ax = Axis(fig[i, j])
            scatter!(ax, ens_gmt, chols[i+shift, j+shift, 1, :], color=:orange, alpha=0.5)
            fit = [chol_coefs[1, i+shift, j+shift, 2].*x .+ chol_coefs[1, i+shift, j+shift, 1] for x in ens_gmt]
            lines!(ax, ens_gmt, fit, color=:blue)
        end
    end
    display(fig)
end



hfile = h5open("data/process/covs_ssp585_100d.hdf5")
covs = zeros(200, 200, 12, 251)
for i in 1:251
    covs[:,:,:,i] = read(hfile, "covs_$i")
end

covhats = zeros(200, 200,12, 251)
for j in 1:12
    for (i,x) in enumerate(ens_gmt)
        U = chol_coefs[j, :, :, 2].*x .+ chol_coefs[j, :, :, 1]
        covhats[:,:,j,i] = U' * U
    end
end



begin
    shift = 190
    fig = Figure(resolution=(1000, 1000))
    for i in 1:5
        for j in 1:5
            ax = Axis(fig[i, j])
            scatter!(ax, ens_gmt, covs[i+shift, j, 8, :], color=:orange, alpha=0.5)

            lines!(ax, ens_gmt, covhats[i+shift, j, 8, :], color=:blue)

            # lines!(ax, ens_gmt, fit, color=:blue)
        end
    end
    display(fig)
end

using_two = true 
second_var ="hurs" # "pr", "huss", "hurs (first variable is always tas)
non_dim = false  
use_metrics = false