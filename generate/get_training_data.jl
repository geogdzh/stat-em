function get_training_data(d::Int)
    scenario = "ssp585"

    ##### assemble training data from projts history
    hfile = h5open("data/$(parent_folder)/projts_historical_$(d)d_$(num_ens_members)ens.hdf5", "r")
    projts_history = read(hfile, "projts")
    ens_gmt_history = read(hfile, "ens_gmt")
    close(hfile)

    hfile = h5open("data/$(parent_folder)/projts_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "r")
    projts_ssp585 = read(hfile, "projts")
    ens_gmt_ssp585 = read(hfile, "ens_gmt")
    close(hfile)

    projts = hcat(projts_history, projts_ssp585)
    ens_gmt = hcat(ens_gmt_history, ens_gmt_ssp585)

    hfile = h5open("data/$(parent_folder)/training_data_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w")
    write(hfile, "projts", projts)
    write(hfile, "ens_gmt", ens_gmt)
    close(hfile)
end