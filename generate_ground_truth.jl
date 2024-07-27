using HDF5
include("./utils/data_util.jl")

#create directories
isdir(pwd() * "/data") ? nothing : mkdir(pwd() * "/data")
isdir(pwd() * "/figs") ? nothing : mkdir(pwd() * "/figs")
isdir(pwd() * "/data/process") ? nothing : mkdir(pwd() * "/data/process")

#specify which variable to add on to temp 
second_var = "huss" # "pr", "huss", "hurs (first variable is always tas)

#get baseline gmt sequences and true 
include("generate/get_gmts.jl")

isdir(pwd() * "/data/ground_truth") ? nothing : mkdir(pwd() * "/data/ground_truth")
include("generate/get_true_var.jl")
