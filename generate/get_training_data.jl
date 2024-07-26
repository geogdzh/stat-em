using HDF5#, ProgressBars
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl")

L1, L2 = 1980, 1032 #for CMIP6
using_two = (ARGS[2] == "true" )
second_var = ARGS[3] # "pr" or "huss"
non_dim = (ARGS[4] == "true" )  
use_metrics = (ARGS[5] == "true" )
if using_two
    if second_var == "pr"
        parent_folder = "temp_precip"
    else
        parent_folder = "temp_huss"
    end
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
d = parse(Int, ARGS[1])
scenario = "ssp585"

##### assemble training data from projts history
hfile = h5open("data/$(parent_folder)/projts_historical_$(d)d_49ens.hdf5", "r")
projts_history = read(hfile, "projts")
ens_gmt_history = read(hfile, "ens_gmt")
close(hfile)

hfile = h5open("data/$(parent_folder)/projts_$(scenario)_$(d)d_50ens.hdf5", "r")
projts_ssp585 = read(hfile, "projts")
ens_gmt_ssp585 = read(hfile, "ens_gmt")
close(hfile)
projts_ssp585 = projts_ssp585[:, :, 1:end .!= 8] #HARDCODED for an error in the 8th historical ens member
ens_gmt_ssp585 = ens_gmt_ssp585[1:end .!= 8, :]

num_ens_members = 49 # number of model runs used to train the emulator
projts = hcat(projts_history, projts_ssp585)
ens_gmt = hcat(ens_gmt_history, ens_gmt_ssp585)


hfile = h5open("data/$(parent_folder)/training_data_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w")
write(hfile, "projts", projts)
write(hfile, "ens_gmt", ens_gmt)
write(hfile, "num_ens_members", num_ens_members)
close(hfile)