using HDF5
include("./utils/data_util.jl")
include("./utils/eof_util.jl")
include("./utils/emulator_util.jl")

#create directories
isdir(pwd() * "/data") ? nothing : mkdir(pwd() * "/data")

#specify which variable to add on to temp 
second_var = "hurs" # "pr", "huss", "hurs (first variable is always tas)

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
L1, L2 = 1980, 1032 #time lengths of historical and future runs
l1, l2 = 165, 86 #num years
M, N = 192, 96
#
ensemble_members = [x for x in 1:30]
deleteat!(ensemble_members, findall(x->x==8,ensemble_members)) #issue in historical tas data
deleteat!(ensemble_members, findall(x->x==3,ensemble_members)) #issue in ssp245 hurs data
num_ens_members = length(ensemble_members)

#get baseline gmt sequences and true 
if !isfile("data/ground_truth/historical_gmts.jld") ## generalize! to make sure they all exist
    include("generate/get_gmts.jl")
end

isdir(pwd() * "/data/ground_truth") ? nothing : mkdir(pwd() * "/data/ground_truth")
include("generate/get_true_var.jl")
