using LinearAlgebra, HDF5, ProgressBars
include("./utils/data_util.jl")
include("./utils/eof_util.jl")
include("./utils/emulator_util.jl")

#specify flags to use
using_two = true 
second_var = "huss" # "pr", "huss", "hurs (first variable is always tas)
non_dim = false  
use_metrics = true
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
isdir(pwd() * "/data/$parent_folder") ? nothing : mkdir(pwd() * "/data/$parent_folder")
isdir(pwd() * "/figs/$parent_folder") ? nothing : mkdir(pwd() * "/figs/$parent_folder")

# specify parameters of data 
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
ensemble_members = [x for x in 1:30]
deleteat!(ensemble_members, findall(x->x==8,ensemble_members)) #issue in historical tas data
deleteat!(ensemble_members, findall(x->x==3,ensemble_members)) #issue in ssp245 hurs data
num_ens_members = length(ensemble_members)
L1, L2 = 1980, 1032 #time lengths of historical and future runs
l1, l2 = 165, 86 #num years
M, N = 192, 96
time_history = [x for x in 1850:2014]
time_future = [x for x in 2015:2100]

#use a sample file to get lat/lon vectors
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.

#generate basis
include("generate/get_basis.jl")

#generate two emulators
for d in [10, 100]   #currently later analysis scripts are hardcoded for just these values of d
    #generate projected timeseries + combine it into training data
    include("generate/get_projts_history.jl")
    include("generate/get_training_data.jl")

    #train emulator
    include("generate/get_emulator.jl")

end

#test the emulators
for param in ["d", "k"]
    if param == "d"
        for d in [10, 100]
            include("generate/get_ens_var.jl")
        end
    elseif param == "k"
        d = 100
        include("generate/get_ens_var.jl")
    end
end
include("generate/get_rmse.jl")

