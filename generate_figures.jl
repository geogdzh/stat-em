using CairoMakie, ProgressBars, HDF5, GeoMakie
include("utils/data_util.jl")
include("utils/eof_util.jl")
include("utils/emulator_util.jl") 

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
isdir(pwd() * "/figs") ? nothing : mkdir(pwd() * "/figs")
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
var_labels = Dict("tas" => "temperature", "huss" => "specific humidity", "hurs" => "relative humidity", "pr" => "precipitation")
scenario_colors = Dict("historical" => :red4, "ssp585" => :red, "ssp245" => :magenta3, "ssp119" => :indigo)
unit_labels = Dict("tas" => "K", "huss" => "kg/kg", "hurs" => "%", "pr" => "mm/day")


######## begin plotting
# include("visualize/figs_data.jl")

# include("visualize/figs_fits.jl")

##
include("visualize/select_sample_locations.jl")

##
include("visualize/figs_samples.jl")
for variable in ["tas", second_var]
    visualize_samples(variable)
end

##
include("visualize/figs_rmse.jl")
for variable in ["tas", second_var]
    generate_rmse_fig(variable)
end

#= implement two potential other plots
adapting from sample locations ==> parallel plots but for a different scenario (validate power of PDFs more) ## and this also feels natural
+ explore_seasonality ==> sample changes in max difference (show advantage of looking at months) ## maybe not this, maybe in an appendix
+ explore_trajectories ==> sample trajectory variablity (show for transition section)=# ## I mean this is definitely necessry