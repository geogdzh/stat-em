using LinearAlgebra, HDF5, ProgressBars
include("./utils/data_util.jl")
include("./utils/eof_util.jl")
include("./utils/emulator_util.jl")

#specify flags to use
hurs_option = "log"

using_two = (ARGS[1] == "true" ) #true 
second_var = ARGS[2] #"hurs" # "pr", "huss", "hurs (first variable is always tas)
non_dim = (ARGS[3] == "true" ) #false  
use_metrics = (ARGS[4] == "true" )  #false
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
isdir(pwd() * "/data/process") ? nothing : mkdir(pwd() * "/data/process")
isdir(pwd() * "/data/$parent_folder") ? nothing : mkdir(pwd() * "/data/$parent_folder")

# specify parameters of data 
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
ensemble_members = [x for x in 1:30] #yeah this needs to carry through somehow
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

println("getting basis")
flush(stdout)
#generate basis
if !isfile("data/$parent_folder/basis_2000d.hdf5")
    include("generate/get_basis.jl")
end

#generate two emulators
include("generate/get_projts_history.jl")
include("generate/get_training_data.jl")
include("generate/get_emulator.jl")

for d in [10, 100] 
    #generate projected timeseries + combine it into training data
    println("working on proj history for emulator with $d modes") 
    flush(stdout)
    get_projts_history(d) # has built in check for historical and ssp585 already exsting - generalize!
    get_training_data(d) #built in check for trainging file
    
    println("training emulator with $d modes")
    flush(stdout)
    #train emulator
    get_emulator(d) #built in checks for emulator already existing
end


#test the emulators
isdir(pwd() * "/data/$parent_folder/ens_vars") ? nothing : mkdir(pwd() * "/data/$parent_folder/ens_vars")
include("generate/get_ens_var.jl")
for param in ["d", "k"]
    if param == "d"
        for d in [10, 100]
            println("calculating ens var for emulator with $d modes")
            flush(stdout)
            run_ens_vars(param, d)
        end
    elseif param == "k"
        d = 100
        run_ens_vars(param, d)
    end
end

println("calculating RMSE")
flush(stdout)
include("generate/get_rmse.jl")
for variable in ["tas", second_var]
    numbers = [10, 100]
    calculate_rmse(numbers, variable, scenarios; for_k=false)
    
    ks = [x for x in 1:2]
    calculate_rmse(ks, variable, scenarios; for_k=true)
end

## generate pattern scaling comparsion
isdir(pwd() * "/data/pattern_scaling") ? nothing : mkdir(pwd() * "/data/pattern_scaling")
include("generate/pattern_scaling.jl")
for variable in ["tas", second_var]
    for measure in ["mean", "std"]
        generate_pattern_scaling(variable, measure)
    end
end