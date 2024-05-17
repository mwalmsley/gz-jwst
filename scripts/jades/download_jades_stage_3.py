import logging

import tqdm
import pandas as pd

from gzjwst.survey_products import get_observations, download_mast_products_simplified
# from astroquery.mast import Observations

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    # PROPOSAL_ID = '1727'
    DATA_DIRECTORY = 'data/jades/stage_3_science_products'

    df = pd.DataFrame()

    # listed on https://archive.stsci.edu/hlsp/jades
    PROPOSAL_IDS = [1180, 1181, 1210, 1286, 1895, 1963, 3215]
    for proposal_id in tqdm.tqdm(PROPOSAL_IDS, unit='proposal'):
        # get all data products
        # for cosmos, this takes about 5 mins, not instant like the others.
        proposal_df = get_observations(proposal_id=proposal_id, calib_level=3)

        # COSMOS Webb has also uploaded all their source catalogs, let's filter those out for now
        # also filter out everything that isn't a miniumum recommended product (mrp)
        # data_products_fits = Observations.filter_products(data_products, extension='fits', mrp_only=True)

        # also omit the auxillary data
        df = pd.concat([df, proposal_df])

    # apply our standard filters
    df = df.query('productGroupDescription == "Minimum Recommended Products"')
    df = df[df['productFilename'].str.endswith('.fits')]
    df = df[df['productType']=='SCIENCE']

    # convert to pandas and save for future reference/debugging
    df.to_csv(DATA_DIRECTORY + '/products_to_download.csv', index=False)  
    assert not any(df['obsID'].duplicated()), "There are duplicate obsIDs in the COSMOS data"

    # TODO set max_files=None to download all
    download_mast_products_simplified(product_uris=df['dataURI'], save_dir=DATA_DIRECTORY, max_files=2, cache=True)  
