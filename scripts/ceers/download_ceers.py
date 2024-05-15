import os
import logging

from gzjwst.survey_products import get_observations, download_mast_products_simplified


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)


    PROVENANCE_NAME = 'ceers'
    DATA_DIRECTORY = os.path.abspath(os.path.dirname(__file__)) + '/../../data/ceers/hlsp'
    # NIRCAM_FILTERS = ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']

    # get all data products associated with CEERS
    data_products = get_observations(provenance_name=PROVENANCE_NAME, query_kwargs={'instrument_name':'NIRCAM/IMAGE'})

    # For pointing 2, there are two sets of observations for F200W and F444W, obtained 1 week apart. 
    # The second set ("b") was affected by significant persistence. 
    # The team provides mosaics of each observation separately (i.e., "2", "2b") as well as combined together ("2all").
    # I'm assuming we only want to keep the 2all images.
    # The easiest way for me to do this was to just note the obsids we want to remove
    additional_pointing2_ids = ['149079958','149074352','149077677','149072294']
    omit_additonal_pointing2_data = [i for i in range(len(data_products)) if data_products['obsID'][i] not in additional_pointing2_ids]
    data_products_fixed_pointing2 = data_products[omit_additonal_pointing2_data]

    # ok, we have 7 filters, 4 pointings, so we should expect 28
    assert len(data_products_fixed_pointing2) == 28, "CEERS should have 28 images, ask Hayley what went wrong"

    # convert to pandas and save for future reference/debugging
    df = data_products_fixed_pointing2.to_pandas()
    assert not any(df['obsID'].duplicated()), "There are duplicate obsIDs in the CEERS data"
    df.to_csv(DATA_DIRECTORY + '/products_to_download.csv', index=False)  

    # TODO set max_files=None to download all
    download_mast_products_simplified(product_uris=df['dataURI'], save_dir=DATA_DIRECTORY, max_files=2, cache=True)  
