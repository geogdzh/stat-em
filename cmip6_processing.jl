#usage: pass in which scenario (number in list) you want to process as an argument

using NCDatasets, HDF5, ProgressBars, DataStructures
include("utils/data_util.jl")

filehead =  "/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/"

scenarios = ["historical", "ssp585", "ssp370", "ssp245", "ssp126", "ssp119"]
scenario = scenarios[parse(Int64, ARGS[1])]
println("processing scenario: $scenario")
variables = ["tas", "pr", "huss"]
variable = variables[parse(Int64, ARGS[2])]
println("processing variables: $variable")
flush(stdout)
sample_year = scenario == "historical" ? 1850 : 2015

sample_filetail = "r$(1)i1p1f1/$(scenario)/$(variable)/250_km/mon/$(sample_year)/CMIP6_MPI-ESM1-2-LR_r$(1)i1p1f1_$(scenario)_$(variable)_250_km_mon_gn_$(sample_year).nc" #for historical tas
sample_file = filehead * sample_filetail
sample_ds = Dataset(sample_file)
latvec = sample_ds["lat"][:]
lonvec = sample_ds["lon"][:]
experiment = sample_ds.attrib["experiment"]
activity_id = sample_ds.attrib["activity_id"]
tracking_id = sample_ds.attrib["tracking_id"]
parent_experiment_id = sample_ds.attrib["parent_experiment_id"]
branch_time_in_child = sample_ds.attrib["branch_time_in_child"]
branch_time_in_parent = sample_ds.attrib["branch_time_in_parent"]
close(sample_ds)
M = length(lonvec)
N = length(latvec)

###
for ens_ind in 1:50
    println("working on ensemble member $ens_ind")
    flush(stdout)
    years = scenario == "historical" ? [x for x in 1850:2014] : [x for x in 2015:2100]
    fulldata = fill(missing, (M, N, length(years)*12))
    fulldata = Array{Union{Missing, Float32}}(fulldata)
    fulltime = Vector{DateTime}(undef, length(years)*12)
    for (i, year) in enumerate(years)
        filetail = "r$(ens_ind)i1p1f1/$(scenario)/$(variable)/250_km/mon/$year/CMIP6_MPI-ESM1-2-LR_r$(ens_ind)i1p1f1_$(scenario)_$(variable)_250_km_mon_gn_$year.nc" 
        file = filehead * filetail
        ds = Dataset(file)
        fulldata[:,:,12*(i-1)+1:12*i] = ds[variable][:,:,:]
        fulltime[12*(i-1)+1:12*i] = ds["time"][:]
        close(ds)
    end

    newnc = "/net/fs06/d3/mgeo/CMIP6/interim/$(scenario)/$(variable)/r$(ens_ind)i1p1f1_$(scenario)_$(variable).nc"
    ds = NCDataset(newnc,"c", attrib = OrderedDict(
        "Conventions"               => "CF-1.7 CMIP-6.2",
        "activity_id"               => activity_id,
        "branch_method"             => "standard",
        "branch_time_in_child"      => branch_time_in_child,
        "branch_time_in_parent"     => branch_time_in_parent,
        "contact"                   => "cmip6-mpi-esm@dkrz.de",
        "creation_date"             => "2019; messily reorganized by mgeo in 2024",
        "data_specs_version"        => "01.00.30",
        "experiment"                => experiment,
        "experiment_id"             => scenario,
        "external_variables"        => "areacella",
        "forcing_index"             => Int32(1),
        "frequency"                 => "mon",
        "further_info_url"          => "https://furtherinfo.es-doc.org/CMIP6.MPI-M.MPI-ESM1-2-LR.$(scenario).none.r$(ens_ind)i1p1f1",
        "grid"                      => "gn",
        "grid_label"                => "gn",
        "history"                   => "2019; CMOR rewrote data to be consistent with CMIP6, CF-1.7 CMIP-6.2 and CF standards.",
        "initialization_index"      => Int32(1),
        "institution"               => "Max Planck Institute for Meteorology, Hamburg 20146, Germany",
        "institution_id"            => "MPI-M",
        "mip_era"                   => "CMIP6",
        "nominal_resolution"        => "250 km",
        "parent_activity_id"        => "CMIP",
        "parent_experiment_id"      => parent_experiment_id,
        "parent_mip_era"            => "CMIP6",
        "parent_source_id"          => "MPI-ESM1-2-LR",
        "parent_time_units"         => "days since 1850-1-1 00:00:00",
        "parent_variant_label"      => "r1i1p1f1",
        "physics_index"             => Int32(1),
        "product"                   => "model-output",
        "project_id"                => "CMIP6",
        "realization_index"         => Int32(ens_ind),
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
        "table_id"                  => "Amon",
        "table_info"                => "Creation Date:(09 May 2019) MD5:e6ef8ececc8f338646ebfb3aeed36bfc",
        "title"                     => "MPI-ESM1-2-LR output prepared for CMIP6",
        "variable_id"               => variable,
        "variant_label"             => "r$(ens_ind)i1p1f1",
        "license"                   => "CMIP6 model data produced by MPI-M is licensed under a Creative Commons Attribution ShareAlike 4.0 International License (https://creativecommons.org/licenses). Consult https://pcmdi.llnl.gov/CMIP6/TermsOfUse for terms of use governing CMIP6 output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further_info_url (recorded as a global attribute in this file) and. The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law.",
        "cmor_version"              => "3.5.0",
        "tracking_id"               => tracking_id,
        "DODS_EXTRA.Unlimited_Dimension" => "time",
    ))

    # Dimensions

    ds.dim["time"] = Inf # unlimited dimension
    ds.dim["bnds"] = 2
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
        "coordinates"               => "height",
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
        "coordinates"               => "height",
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
        "coordinates"               => "height",
    ))

    if variable == "tas"
        ncheight = defVar(ds,"height", Float64, (), attrib = OrderedDict(
            "_FillValue"                => NaN,
            "units"                     => "m",
            "axis"                      => "Z",
            "positive"                  => "up",
            "long_name"                 => "height",
            "standard_name"             => "height",
        ))

        nctas = defVar(ds,"tas", Float32, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float32(1.0e20),
            "standard_name"             => "air_temperature",
            "long_name"                 => "Near-Surface Air Temperature",
            "comment"                   => "near-surface (usually, 2 meter) air temperature",
            "units"                     => "K",
            "cell_methods"              => "area: time: mean",
            "cell_measures"             => "area: areacella",
            "history"                   => "2019-09-12T09:26:24Z altered by CMOR: Treated scalar dimension: 'height'. 2019-09-12T09:26:24Z altered by CMOR: replaced missing value flag (-9e+33) and corresponding data with standard missing value (1e+20). 2019-09-12T09:26:24Z altered by CMOR: Inverted axis: lat.",
            "_ChunkSizes"               => Int32[1, 96, 192],
            "coordinates"               => "height",
            "missing_value"             => Float32(1.0e20),
        ))

        #define the data
        # ncheight[:] = ...
        nctas[:] = fulldata

    elseif variable == "pr"
        ncpr = defVar(ds,"pr", Float32, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float32(1.0e20),
            "standard_name"             => "precipitation_flux",
            "long_name"                 => "Precipitation",
            "comment"                   => "includes both liquid and solid phases",
            "units"                     => "kg m-2 s-1",
            "original_name"             => "pr",
            "cell_methods"              => "area: time: mean",
            "cell_measures"             => "area: areacella",
            "history"                   => "2019-09-11T14:13:17Z altered by CMOR: replaced missing value flag (-9e+33) and corresponding data with standard missing value (1e+20). 2019-09-11T14:13:18Z altered by CMOR: Inverted axis: lat.",
            "_ChunkSizes"               => Int32[1, 96, 192],
            "missing_value"             => Float32(1.0e20),
        ))

        #define the data
        ncpr[:] = fulldata

    elseif variable == "huss"
        nchuss = defVar(ds,"huss", fulldata, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float32(1.0e20),
            "standard_name"             => "specific_humidity",
            "long_name"                 => "Near-Surface Specific Humidity",
            "comment"                   => "Near-surface (usually, 2 meter) specific humidity.",
            "units"                     => "1",
            "cell_methods"              => "area: time: mean",
            "cell_measures"             => "area: areacella",
            "history"                   => "2019-09-11T14:13:17Z altered by CMOR: Treated scalar dimension: 'height'. 2019-09-11T14:13:17Z altered by CMOR: replaced missing value flag (-9e+33) and corresponding data with standard missing value (1e+20). 2019-09-11T14:13:18Z altered by CMOR: Inverted axis: lat.",
            "_ChunkSizes"               => Int32[1, 96, 192],
            "coordinates"               => "height",
            "missing_value"             => Float32(1.0e20),
        ))

        # nchuss[:,:,:] = fulldata
    end

    # Define variables

    nctime[:] = fulltime
    # nctime_bnds[:] = ...
    nclat[:] = latvec
    # nclat_bnds[:] = ...
    nclon[:] = lonvec
    # nclon_bnds[:] = ...
    

    close(ds)

end
