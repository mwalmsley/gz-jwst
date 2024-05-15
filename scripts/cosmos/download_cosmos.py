from gzjwst.survey_products import get_observations, download_mast_products
from astroquery.mast import Observations

PROPOSAL_ID = '1727'
DATA_DIRECTORY = 'data/cosmos/hlsp'

# get all data products
data_products = get_observations(proposal_id=PROPOSAL_ID, query_kwargs={'instrument_name':'NIRCAM/IMAGE'})

# COSMOS Webb has also uploaded all their source catalogs, let's filter those out for now
# also filter out everything that isn't a miniumum recommended product (mrp)
data_products_fits = Observations.filter_products(data_products, extension='fits', mrp_only=True)

# also omit the auxillary data
data_products_fits = data_products_fits[data_products_fits['productType']=='SCIENCE']


# test call - downloads just a script
download_mast_products(data_products_fits, save_dir=DATA_DIRECTORY,  curl_flag=True)

# real script
# download_mast_products(data_products_fits, save_dir=DATA_DIRECTORY)