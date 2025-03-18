import os
import logging

import pandas as pd

from gzjwst.survey_products import get_observations, download_mast_products_simplified


def make_ceers_product_list():

    # get all data products associated with CEERS
    data_products = get_observations(provenance_name=PROVENANCE_NAME)

    # For pointing 2, there are two sets of observations for F200W and F444W, obtained 1 week apart. 
    # The second set ("b") was affected by significant persistence. 
    # The team provides mosaics of each observation separately (i.e., "2", "2b") as well as combined together ("2all").
    # I'm assuming we only want to keep the 2all images.
    # The easiest way to do this was to just note the obsids we want to remove
    additional_pointing2_ids = ['149079958','149074352','149077677','149072294']
    df = data_products[~data_products['obsID'].isin(additional_pointing2_ids)]

    # ok, we have 7 filters, 4 pointings, so we should expect 28
    assert len(df) == 28, "CEERS should have 28 images, ask Hayley what went wrong"

    #  save for future reference/debugging
    assert not any(df['obsID'].duplicated()), "There are duplicate obsIDs in the CEERS data"
    return df


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)


    # this uses the CEERS HLSP products
    PROVENANCE_NAME = 'ceers'
    DATA_DIRECTORY = 'data/ceers/hlsp'


    df = make_ceers_product_list()
    df.to_csv(DATA_DIRECTORY + '/products_to_download.csv', index=False)  

    df = pd.read_csv('data/ceers/hlsp/products_to_download.csv')  # everything
    # df = pd.read_csv('data/ceers/hlsp/products_to_download_nircam3_only.csv')  # temporary, for SEP testing

    # TODO set max_files=None to download all
    download_mast_products_simplified(product_uris=df['dataURI'], save_dir=DATA_DIRECTORY, max_files=10, cache=True)  
