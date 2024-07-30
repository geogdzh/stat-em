variables = ["tas", second_var]  #could also call this from the outer script

for variable in variables
    println("working on $(variable)")
    flush(stdout)
    for scenario in scenarios   
        println("working on $(scenario)")
        flush(stdout)
        if isfile("data/ground_truth/vars_$(variable)_$(scenario)_$(num_ens_members)ens.hdf5")
            println("vars for $(variable) and $(scenario) already exist")
            flush(stdout)
            continue
        end

        L = scenario == "historical" ? L1 : L2

        all_gmts = zeros((num_ens_members, Int(L/12))) 
        all_data = zeros((M,N, L, num_ens_members))
        for (i, n) in enumerate(ensemble_members)
            println("working on ensemble member $(n)")
            flush(stdout)
            file = file_head*"$(scenario)/$(variable)/r$(n)i1p1f1_$(scenario)_$(variable).nc"
            ts = ncData(file, variable) 
            all_gmts[i,:] = get_gmt_list(ts)
            all_data[:,:,:,i] = apply_transform(ts.data[:,:,:], variable; hurs_option=hurs_option)
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