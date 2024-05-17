import logging

from gzjwst.survey_products import get_observations, download_mast_products_simplified
from astroquery.mast import Observations

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    PROPOSAL_ID = '1727'
    DATA_DIRECTORY = 'data/cosmos/stage_3_science_products'

    # get all data products
    # for cosmos, this takes about 5 mins, not instant like the others.
    df = get_observations(proposal_id=PROPOSAL_ID, calib_level=3)

    # apply our standard filters
    df = df.query('productGroupDescription == "Minimum Recommended Products"')
    df = df[df['productFilename'].str.endswith('.fits')]
    df = df[df['productType']=='SCIENCE']

    # convert to pandas and save for future reference/debugging
    df.to_csv(DATA_DIRECTORY + '/products_to_download.csv', index=False)  
    assert not any(df['obsID'].duplicated()), "There are duplicate obsIDs in the COSMOS data"

    # TODO set max_files=None to download all
    download_mast_products_simplified(product_uris=df['dataURI'], save_dir=DATA_DIRECTORY, max_files=2, cache=True)  
