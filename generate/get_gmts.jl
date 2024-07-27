for scenario in scenarios
    ens_gmt = scenario == "historical" ? zeros((num_ens_members, Int((L1)/12))) : zeros((num_ens_members, Int((L2)/12))) #GMT in any case

    errors = []
    for i in 1:num_ens_members
        try
            println("working on ensemble member $(i)")
            flush(stdout)
            files = file_head*"$(scenario)/tas/r$(i)i1p1f1_$(scenario)_tas.nc"
            tmps = ncData(files, "tas")
            ens_gmt[i, :] = get_gmt_list(tmps) 
        catch
            println("missing values in ensemble member $(i)")
            push!(errors, i)
            flush(stdout)
        end
    end

    hfile = h5open("data/$(scenario)_gmts_$(num_ens_members)ens.hdf5", "w")
    write(hfile, "ens_gmt", ens_gmt)    
    write(hfile, "num_ens_members", num_ens_members)
    close(hfile)
end