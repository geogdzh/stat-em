
# get the EOF coefficients
hfile = h5open("data/$parent_folder/training_data_ssp585_$(d)d_$(num_ens_members)ens.hdf5", "r")
ens_gmt = read(hfile, "ens_gmt")
ens_projts = read(hfile, "projts")
close(hfile)

