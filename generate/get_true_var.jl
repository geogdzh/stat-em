variables = ["tas", second_var]

for variable in variables
    println("working on $(variable)")
    flush(stdout)
    for scenario in scenarios   
        println("working on $(scenario)")

        L = scenario == "historical" ? L1 : L2

        all_gmts = zeros((num_ens_members, Int(L/12))) 
        all_data = zeros((M,N, L, num_ens_members))
        for i in 1:num_ens_members
            println("working on ensemble member $(i)")
            flush(stdout)
            file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
            file = file_head*"$(scenario)/$(variable)/r$(i)i1p1f1_$(scenario)_$(variable).nc"
            ts = ncData(file, variable) 
            all_gmts[i,:] = get_gmt_list(ts)
            all_data[:,:,:,i] = apply_transform(ts.data[:,:,:], variable)
        end

        true_var = var(all_data, dims=4)[:,:,:,1]
        true_ens_gmt = mean(all_gmts, dims=1)[:]
        true_ens_mean = mean(all_data, dims=4)[:,:,:,1]

        hfile = h5open("data/ground_truth/vars_$(variable)_$(scenario)_$(num_ens_members)ens.hdf5", "w")
        write(hfile, "true_var", true_var)
        write(hfile, "true_ens_gmt", true_ens_gmt)
        write(hfile, "true_ens_mean", true_ens_mean)
        close(hfile)
    end
end