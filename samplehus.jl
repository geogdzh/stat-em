using NCDatasets, DataStructures
ds = NCDataset("filename.nc","c", attrib = OrderedDict(
    "Conventions"               => "CF-1.7 CMIP-6.2",
    "activity_id"               => "CMIP",
    "branch_method"             => "standard",
    "branch_time_in_child"      => 0.0,
    "branch_time_in_parent"     => 0.0,
    "contact"                   => "cmip6-mpi-esm@dkrz.de",
    "creation_date"             => "2019-09-11T14:13:25Z",
    "data_specs_version"        => "01.00.30",
    "experiment"                => "all-forcing simulation of the recent past",
    "experiment_id"             => "historical",
    "external_variables"        => "areacella",
    "forcing_index"             => Int32(1),
    "frequency"                 => "mon",
    "further_info_url"          => "https://furtherinfo.es-doc.org/CMIP6.MPI-M.MPI-ESM1-2-LR.historical.none.r1i1p1f1",
    "grid"                      => "gn",
    "grid_label"                => "gn",
    "history"                   => "2019-09-11T14:13:25Z ; CMOR rewrote data to be consistent with CMIP6, CF-1.7 CMIP-6.2 and CF standards.",
    "initialization_index"      => Int32(1),
    "institution"               => "Max Planck Institute for Meteorology, Hamburg 20146, Germany",
    "institution_id"            => "MPI-M",
    "mip_era"                   => "CMIP6",
    "nominal_resolution"        => "250 km",
    "parent_activity_id"        => "CMIP",
    "parent_experiment_id"      => "piControl",
    "parent_mip_era"            => "CMIP6",
    "parent_source_id"          => "MPI-ESM1-2-LR",
    "parent_time_units"         => "days since 1850-1-1 00:00:00",
    "parent_variant_label"      => "r1i1p1f1",
    "physics_index"             => Int32(1),
    "product"                   => "model-output",
    "project_id"                => "CMIP6",
    "realization_index"         => Int32(1),
    "realm"                     => "atmos",
    "references"                => "MPI-ESM: Mauritsen, T. et al. (2019), Developments in the MPI‐M Earth System Model version 1.2 (MPI‐ESM1.2) and Its Response to Increasing CO2, J. Adv. Model. Earth Syst.,11, 998-1038, doi:10.1029/2018MS001400,
Mueller, W.A. et al. (2018): A high‐resolution version of the Max Planck Institute Earth System Model MPI‐ESM1.2‐HR. J. Adv. Model. EarthSyst.,10,1383–1413, doi:10.1029/2017MS001217",
    "source"                    => "MPI-ESM1.2-LR (2017): 
aerosol: none, prescribed MACv2-SP
atmos: ECHAM6.3 (spectral T63; 192 x 96 longitude/latitude; 47 levels; top level 0.01 hPa)
atmosChem: none
land: JSBACH3.20
landIce: none/prescribed
ocean: MPIOM1.63 (bipolar GR1.5, approximately 1.5deg; 256 x 220 longitude/latitude; 40 levels; top grid cell 0-12 m)
ocnBgchem: HAMOCC6
seaIce: unnamed (thermodynamic (Semtner zero-layer) dynamic (Hibler 79) sea ice model)",
    "source_id"                 => "MPI-ESM1-2-LR",
    "source_type"               => "AOGCM",
    "sub_experiment"            => "none",
    "sub_experiment_id"         => "none",
    "table_id"                  => "Emon",
    "table_info"                => "Creation Date:(09 May 2019) MD5:e6ef8ececc8f338646ebfb3aeed36bfc",
    "title"                     => "MPI-ESM1-2-LR output prepared for CMIP6",
    "variable_id"               => "hus",
    "variant_label"             => "r1i1p1f1",
    "license"                   => "CMIP6 model data produced by MPI-M is licensed under a Creative Commons Attribution ShareAlike 4.0 International License (https://creativecommons.org/licenses). Consult https://pcmdi.llnl.gov/CMIP6/TermsOfUse for terms of use governing CMIP6 output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further_info_url (recorded as a global attribute in this file) and. The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law.",
    "cmor_version"              => "3.5.0",
    "tracking_id"               => "hdl:21.14100/4ce1f2f8-94b0-4fbb-9db7-473041b61281",
    "DODS_EXTRA.Unlimited_Dimension" => "time",
))

# Dimensions

ds.dim["time"] = Inf # unlimited dimension
ds.dim["bnds"] = 2
ds.dim["plev"] = 28
ds.dim["lat"] = 96
ds.dim["lon"] = 192

# Declare variables

nctime = defVar(ds,"time", Float64, ("time",), attrib = OrderedDict(
    "_FillValue"                => NaN,
    "bounds"                    => "time_bnds",
    "axis"                      => "T",
    "long_name"                 => "time",
    "standard_name"             => "time",
    "_ChunkSizes"               => Int32(1),
    "units"                     => "days since 1850-01-01",
    "calendar"                  => "proleptic_gregorian",
))

nctime_bnds = defVar(ds,"time_bnds", Float64, ("bnds", "time"), attrib = OrderedDict(
    "_FillValue"                => NaN,
    "_ChunkSizes"               => Int32[1, 2],
))

ncplev = defVar(ds,"plev", Float64, ("plev",), attrib = OrderedDict(
    "_FillValue"                => NaN,
    "units"                     => "Pa",
    "axis"                      => "Z",
    "positive"                  => "down",
    "long_name"                 => "pressure",
    "standard_name"             => "air_pressure",
))

nclat = defVar(ds,"lat", Float64, ("lat",), attrib = OrderedDict(
    "_FillValue"                => NaN,
    "bounds"                    => "lat_bnds",
    "units"                     => "degrees_north",
    "axis"                      => "Y",
    "long_name"                 => "Latitude",
    "standard_name"             => "latitude",
))

nclat_bnds = defVar(ds,"lat_bnds", Float64, ("bnds", "lat"), attrib = OrderedDict(
    "_FillValue"                => NaN,
    "_ChunkSizes"               => Int32[96, 2],
))

nclon = defVar(ds,"lon", Float64, ("lon",), attrib = OrderedDict(
    "_FillValue"                => NaN,
    "bounds"                    => "lon_bnds",
    "units"                     => "degrees_east",
    "axis"                      => "X",
    "long_name"                 => "Longitude",
    "standard_name"             => "longitude",
))

nclon_bnds = defVar(ds,"lon_bnds", Float64, ("bnds", "lon"), attrib = OrderedDict(
    "_FillValue"                => NaN,
    "_ChunkSizes"               => Int32[192, 2],
))

nchus = defVar(ds,"hus", Float32, ("lon", "lat", "plev", "time"), attrib = OrderedDict(
    "_FillValue"                => Float32(1.0e20),
    "standard_name"             => "specific_humidity",
    "long_name"                 => "Specific Humidity",
    "comment"                   => "Specific humidity is the mass fraction of water vapor in (moist) air.",
    "units"                     => "1",
    "cell_methods"              => "time: mean",
    "cell_measures"             => "area: areacella",
    "history"                   => "2019-09-11T14:13:25Z altered by CMOR: Reordered dimensions, original order: time lat lon plev. 2019-09-11T14:13:25Z altered by CMOR: replaced missing value flag (-9e+33) and corresponding data with standard missing value (1e+20). 2019-09-11T14:13:25Z altered by CMOR: Inverted axis: lat.",
    "_ChunkSizes"               => Int32[1, 28, 96, 192],
    "missing_value"             => Float32(1.0e20),
))


# Define variables

# nctime[:] = ...
# nctime_bnds[:,:] = ...
# ncplev[:] = ...
# nclat[:] = ...
# nclat_bnds[:,:] = ...
# nclon[:] = ...
# nclon_bnds[:,:] = ...
# nchus[:,:,:,:] = ...

close(ds)
