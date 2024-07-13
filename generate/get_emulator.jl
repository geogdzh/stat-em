using HDF5#, ProgressBars
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl")

L1, L2 = 1980, 1032 #for CMIP6
scenario = "ssp585"
# scenario = "ssp119"
offload = false
using_precip = true 
non_dim = false  
use_metrics = false
if using_precip
    parent_folder = "temp_precip"
else
    parent_folder = "temp"
end
if non_dim
    parent_folder = "nondim"
end
if use_metrics && using_precip
    parent_folder = "metrics"
elseif use_metrics && !using_precip
    parent_folder = "temp_metrics"
end


######### load in basis
hfile = using_precip ? h5open("data/$(parent_folder)/temp_precip_basis_2000d.hdf5", "r") : h5open("data/$(parent_folder)/temp_basis_2000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
if non_dim ##none of this is actually used here
    temp_factor = read(hfile, "temp_factor")
    pr_factor = read(hfile, "pr_factor")
end
if use_metrics
    metric = read(hfile, "metric")
end
close(hfile)
d = parse(Int, ARGS[1])
basis = basis[:, 1:d]

########## load in training data
hfile = using_precip ?  h5open("data/$(parent_folder)/training_data_withpr_$(scenario)_$(d)d_49ens.hdf5", "r") : h5open("data/$(parent_folder)/training_data_$(scenario)_$(d)d_49ens.hdf5", "r")
ens_projts = read(hfile, "projts")[1:d, :, :]
ens_gmt = read(hfile, "ens_gmt")
num_ens_members = read(hfile, "num_ens_members")
close(hfile)

########## get the emulator itself
println("working on emulator for $(scenario) with $(d) dimensions")
flush(stdout)
ens_gmt = mean(ens_gmt, dims=1)
mean_coefs_1 = get_mean_coefs(ens_projts, ens_gmt, degree=1)
mean_coefs_2 = get_mean_coefs(ens_projts, ens_gmt, degree=2)
mean_coefs_3 = get_mean_coefs(ens_projts, ens_gmt, degree=3)
println("getting the covariances")
flush(stdout)
gmt_cov(ens_projts, ens_gmt, "$(scenario)_$(d)d$(non_dim ? "_nondim" : "")") #saves out covs
println("getting the cholesky decomposition and fitting")
flush(stdout)
chol_coefs = get_chol_coefs(ens_gmt, "$(scenario)_$(d)d$(non_dim ? "_nondim" : "")"; offload=offload) # this is NOT IMPLEMENTED YET (offload=true impossible)


hfile = using_precip ? h5open("data/$(parent_folder)/gaussian_emulator_withpr_$(scenario)_$(d)d.hdf5", "w") : h5open("data/$(parent_folder)/gaussian_emulator_$(scenario)_$(d)d.hdf5", "w")
write(hfile, "mean_coefs_1", mean_coefs_1)
write(hfile, "mean_coefs_2", mean_coefs_2)
write(hfile, "mean_coefs_3", mean_coefs_3)
write(hfile, "chol_coefs", chol_coefs)
write(hfile, "basis", basis)
if non_dim
    write(hfile, "temp_factor", temp_factor)
    write(hfile, "pr_factor", pr_factor)
end
if use_metrics
    write(hfile, "metric", metric)
end
write(hfile, "num_ens_members", num_ens_members)
write(hfile, "ens_gmt", ens_gmt)
close(hfile)


# ###

# begin
#     fig = Figure(resolution=(1500, 1000))
#     ax = Axis(fig[1,1])
#     for i in 1:49
#         scatter!(ax, ens_gmt[:], -1 .* ens_projts[20,1:12:end,i], color=:black, alpha=0.3)
#         #[time_history..., time_future...]
#     end
#     means = mean(ens_projts, dims=3)[:,:]
#     scatter!(ax, ens_gmt[:], -1 .* means[20, 1:12:end], color=:red, alpha=1)
#     # plot!(ax, ens_gmt[:], ens_gmt[:], color=:red, alpha=0.5)
#     display(fig)
# end
# mean(month_to_year_avg(means[1,:]) .+ ens_gmt[:])