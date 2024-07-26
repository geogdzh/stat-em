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

############## load a basis
hfile = h5open("data/$(parent_folder)/basis_2000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
if non_dim
    temp_factor = read(hfile, "temp_factor")
    if using_two
        pr_factor = read(hfile, "pr_factor")
    end
end
if use_metrics
    metric = read(hfile, "metric")
end
close(hfile)
d = parse(Int, ARGS[1])
basis = basis[:, 1:d]

scenarios = ["historical", "ssp585", "ssp245", "ssp119"]

############### 
for scenario in scenarios
    println("working on $(scenario)")
    flush(stdout)
    num_ens_members = 50 # number of model runs used to train the emulator
    projts = scenario == "historical" ? zeros((d, (L1), num_ens_members)) : zeros((d, (L2), num_ens_members))
    ens_gmt = scenario == "historical" ? zeros((num_ens_members, Int((L1)/12))) : zeros((num_ens_members, Int((L2)/12)))

    file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
    errors = []
    for i in 1:num_ens_members
        try
            println("working on ensemble member $(i)")
            flush(stdout)
            files = file_head*"$(scenario)/tas/r$(i)i1p1f1_$(scenario)_tas.nc"
            files_pr = file_head*"$(scenario)/$(second_var)/r$(i)i1p1f1_$(scenario)_$(second_var).nc"
            tmps = ncData(files, "tas")
            prs = using_two ? ncData(files_pr, second_var) : nothing
            if using_two
                data1 = use_metrics ? sqrt.(metric) .* tmps.data : tmps.data
                data1 = non_dim ? reshape_data(data1) ./ temp_factor : reshape_data(data1)

                prdata = second_var == "pr" ? log.(prs.data .* 86400) : prs.data #only applies it if using precip
                data2 = use_metrics ? sqrt.(metric) .* prdata : prdata
                data2 = non_dim ? reshape_data(data2) ./ pr_factor : reshape_data(data2)
                data = vcat(data1, data2)
                projts[:, :, i] =  project_timeseries(data, basis, reshaped=true)
            else
                projts[:, :, i] =  non_dim ? project_timeseries(tmps.data, basis) ./ temp_factor : project_timeseries(tmps.data, basis)
                if use_metrics
                    throw("metrics not implemented for single variable")
                end
            end
            ens_gmt[i, :] = get_gmt_list(tmps) 
        catch
            println("missing values in ensemble member $(i)")
            push!(errors, i)
            flush(stdout)
        end
    end

    if length(errors) > 0
        projts = projts[:,:,1:end .!= errors[:]] # HARDCODED FOR ONLY ONE ERROR - need to generalize. issue was scoping in script run # adn for some reason error when no errors?
        ens_gmt = ens_gmt[1:end .!= errors[:], :]
    end
    num_ens_members = size(ens_gmt)[1]

    hfile = h5open("data/$(parent_folder)/projts_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w")
    write(hfile, "projts", projts)
    write(hfile, "ens_gmt", ens_gmt)
    write(hfile, "num_ens_members", num_ens_members)
    if non_dim
        write(hfile, "temp_factor", temp_factor)
        if using_two
            write(hfile, "pr_factor", pr_factor)
        end
    end
    if use_metrics
        write(hfile, "metric", metric)
    end
    close(hfile)
end