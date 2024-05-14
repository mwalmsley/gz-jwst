from gzjwst.survey_products import get_observations_by_provenance_name, download_mast_products

PROVENANCE_NAME = 'ceers'
CEERS_DATA_DIRECTORY = '../../data/ceers/hlsp'
CEERS_NIRCAM_FILTERS = ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']

# get all data products associated with CEERS
data_products = get_observations_by_provenance_name(PROVENANCE_NAME)

# only select for those that are NIRCAM images
nircam_filter = [i for i in range(len(data_products)) if data_products[i]['filters'] in CEERS_NIRCAM_FILTERS]
data_products_nircam = data_products[nircam_filter]


# For pointing 2, there are two sets of observations for F200W and F444W, obtained 1 week apart. 
# The second set ("b") was affected by significant persistence. 
# The team provides mosaics of each observation separately (i.e., "2", "2b") as well as combined together ("2all").
# I'm assuming we only want to keep the 2all images.
# The easiest way for me to do this was to just note the obsids we want to remove
additional_pointing2_ids = ['149079958','149074352','149077677','149072294']
omit_additonal_pointing2_data = [i for i in range(len(data_products_nircam)) if data_products_nircam['obsID'][i] not in additional_pointing2_ids]
data_products_nircam_fixed_pointing2 = data_products_nircam[omit_additonal_pointing2_data]

# ok, we have 7 filters, 4 pointings, so we should expect 28
assert len(data_products_nircam_fixed_pointing2) == 28, "CEERS should have 28 images, ask Hayley what went wrong"

# test call - downloads just a script
download_mast_products(data_products_nircam_fixed_pointing2, save_dir=CEERS_DATA_DIRECTORY, max_files=28, chunk_size=5, curl_flag=True)

# real script
# download_mast_products(data_products_nircam_fixed_pointing2, save_dir=CEERS_DATA_DIRECTORY, max_files=28, chunk_size=5)