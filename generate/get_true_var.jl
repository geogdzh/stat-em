using ProgressBars, HDF5
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl") 

M, N = 192, 96
L1, L2 = 1980, 1032

# scenarios = ["ssp585", "ssp370", "ssp245", "ssp126", "ssp119", "historical"]
scenarios = ["ssp245", "ssp119", "historical",]
global num_ens_members = 50
variables = ["huss"]#["tas", "pr"]

parent_folder = "ground_truth"

for variable in variables
    println("working on $(variable)")
    flush(stdout)
    for scenario in scenarios   
        println("working on $(scenario)")

        L = scenario == "historical" ? L1 : L2

        all_gmts = zeros((num_ens_members, Int(L/12))) 
        all_data = zeros((M,N, L, num_ens_members))
        for i in 1:num_ens_members
            if i == 8 && scenario == "historical"
                println("skipping missing data in ens member 8")
            else
                println("working on ensemble member $(i)")
                flush(stdout)
                file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
                file = file_head*"$(scenario)/$(variable)/r$(i)i1p1f1_$(scenario)_$(variable).nc"
                ts = ncData(file, variable) 
                all_gmts[i,:] = get_gmt_list(ts)
                if variable == "pr"
                    all_data[:,:,:,i] = log.(ts.data[:,:,:] .* 86400) #convert to mm/day and use log
                else
                    all_data[:,:,:,i] = ts.data[:,:,:]
                end
            end
        end
        if scenario == "historical"
            all_data = all_data[:,:,:,1:end .!= 8] # HARDCODED FOR ONLY ONE ERROR - need to generalize. issue was scoping in script run
            all_gmts = all_gmts[1:end .!= 8, :]
            global num_ens_members = 49
        else
            global num_ens_members = 50
        end
        true_var = var(all_data, dims=4)[:,:,:,1]
        true_ens_gmt = mean(all_gmts, dims=1)[:]
        true_ens_mean = mean(all_data, dims=4)[:,:,:,1]

        hfile = h5open("data/$(parent_folder)/vars_$(variable)_$(scenario)_$(num_ens_members)ens.hdf5", "w")
        write(hfile, "true_var", true_var)
        write(hfile, "num_ens_members", num_ens_members)
        write(hfile, "true_ens_gmt", true_ens_gmt)
        write(hfile, "true_ens_mean", true_ens_mean)
        close(hfile)

    end
end