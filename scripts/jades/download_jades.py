import logging

from gzjwst.survey_products import get_observations, download_mast_products_simplified

def check_goods_field(obs_id: str):
    """
    Check if the obsID is from GOODS-S or GOODS-N
    """
    if 'goods-s' in obs_id:
        return 'goods_s'
    elif 'goods-n' in obs_id:
        return 'goods_n'
    else:
        return None

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO) 

    PROVENANCE_NAME = 'jades'
    DATA_DIRECTORY = 'data/jades/hlsp'
    # NIRCAM_FILTERS = ['F090W', 'F115W', 'F150W', 'F182M', 'F200W', 'F210M', 'F277W', 'F335M', 'F356W', 'F410M', 'F444W']

    # get all data products associated with CEERS
    data_products = get_observations(provenance_name=PROVENANCE_NAME, query_kwargs={'instrument_name':'NIRCAM/IMAGE'})

    # convert to pandas and save for future reference/debugging
    df = data_products.to_pandas()
    df['field'] = df['obs_id'].apply(check_goods_field)
    assert not any(df['obs_id'].duplicated()), "There are duplicate obsIDs in the JADES data"
    df.to_csv(DATA_DIRECTORY + '/products_to_download.csv', index=False)

    # now split into goods-s and goods-n data
    # goods_s = [i for i in range(len(data_products)) if 'goods-s' in data_products['obs_id'][i]]
    # goods_n = [i for i in range(len(data_products)) if 'goods-n' in data_products['obs_id'][i]]

    # data_products_goods_s = data_products[goods_s]
    # data_products_goods_n = data_products[goods_n]


    download_mast_products_simplified(df['dataURI'], save_dir=DATA_DIRECTORY, max_files=2, cache=True)

    # I'm not sure we want to put these in different directories based on field, but if so:
    # download_mast_products_simplified(df.query('field == "goods_s')['dataURI'], save_dir=DATA_DIRECTORY+'goods_s/', max_files=2, cache=True)
    # download_mast_products_simplified(df.query('field == "goods_n')['dataURI'], save_dir=DATA_DIRECTORY+'goods_n/', max_files=2, cache=True)
