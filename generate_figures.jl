using CairoMakie, ProgressBars, HDF5, GeoMakie
include("utils/data_util.jl")
include("utils/eof_util.jl")
include("utils/emulator_util.jl") 

#specify flags to use
using_two = true 
second_var = "hurs" # "pr", "huss", "hurs (first variable is always tas)
non_dim = false  
use_metrics = false
if using_two
    parent_folder = "temp_$second_var"
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
end # this is not an exhaustive list of possible combinations, of course... 
isdir(pwd() * "/figs/$parent_folder") ? nothing : mkdir(pwd() * "/figs/$parent_folder")


# get sample lat/lon vectors
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.

#set more parameters
scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
time_history = [x for x in 1850:2014]
time_future = [x for x in 2015:2100]
L1, L2 = 1980, 1032 
l1, l2 = 165, 86
M, N = 192, 96

# read off number of ens members #FIX THIS so it's automated to match what was used in the emulator
ensemble_members = [x for x in 1:30]
deleteat!(ensemble_members, findall(x->x==8,ensemble_members)) #issue in historical tas data
deleteat!(ensemble_members, findall(x->x==3,ensemble_members)) #issue in ssp245 hurs data
num_ens_members = length(ensemble_members)


##
include("visualize/select_sample_locations.jl")

##
# for variable in ["tas", second_var]
#     include("visualize/figs_samples.jl")
# end
