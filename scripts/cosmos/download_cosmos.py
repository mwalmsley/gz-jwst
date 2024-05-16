import logging

from gzjwst.survey_products import get_observations, download_mast_products_simplified
from astroquery.mast import Observations

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    PROPOSAL_ID = '1727'
    DATA_DIRECTORY = 'data/cosmos/stage_3_science_products'

    # get all data products
    data_products = get_observations(proposal_id=PROPOSAL_ID, query_kwargs={'instrument_name':'NIRCAM/IMAGE'})

    # COSMOS Webb has also uploaded all their source catalogs, let's filter those out for now
    # also filter out everything that isn't a miniumum recommended product (mrp)
    data_products_fits = Observations.filter_products(data_products, extension='fits', mrp_only=True)

    # also omit the auxillary data
    data_products_fits = data_products_fits[data_products_fits['productType']=='SCIENCE']

    # convert to pandas and save for future reference/debugging
    df = data_products_fits.to_pandas()
    assert not any(df['obsID'].duplicated()), "There are duplicate obsIDs in the COSMOS data"
    df.to_csv(DATA_DIRECTORY + '/products_to_download.csv', index=False)  

    # TODO set max_files=None to download all
    download_mast_products_simplified(product_uris=df['dataURI'], save_dir=DATA_DIRECTORY, max_files=2, cache=True)  
