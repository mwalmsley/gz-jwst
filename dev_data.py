import pandas as pd

def make_new_dev_data_product_list():
    # run a test download for each survey first to create these csvs

    ceers = pd.read_csv('data/ceers/hlsp/products_to_download.csv')
    # jades = pd.read_csv('data/jades/hlsp/products_to_download.csv')
    jades = pd.read_csv('data/jades/stage_3_science_products/products_to_download.csv')
    cosmos = pd.read_csv('data/cosmos/stage_3_science_products/products_to_download.csv')

    ceers['survey'] = 'ceers'
    jades['survey'] = 'jades'
    cosmos['survey'] = 'cosmos'

    # ceers is per pointing
    # I picked a random pointing (6)
    ceers_examples = ceers[ceers['obs_id'].str.contains('nircam_nircam6')]

    # jades is not per pointing, it is fully reduced
    # jades_examples = jades[jades['obs_id'].str.contains('goods-s')]  # HLSP mosaic
    jades_examples = jades[jades['obs_id'].str.contains('o006_t008')]

    # cosmos is per pointing
    # I picked a random pointing (140)
    cosmos_examples = cosmos[cosmos['obs_id'].str.contains('o140_t104')]
    # print(cosmos_examples['filters'].values)

    # only wide filters?
    ceers_examples = ceers_examples[ceers_examples['filters'].str.contains('W')]
    jades_examples = jades_examples[jades_examples['filters'].str.contains('W')]

    df = pd.concat([ceers_examples, jades_examples, cosmos_examples])

    print(df['size'].sum() / 1e9) # GB
    print(df.groupby('survey').agg({'size': 'sum'}) / 1e9) # GB
    """
    size
    survey          
    ceers   8.101352
    cosmos  3.640905
    jades   6.929211
    """

    df.to_csv('data//dev_data/dev_data.csv', index=False)  # important columns are survey and filters

if __name__ == '__main__':

    # make_new_dev_data_product_list()

    df = pd.read_csv('data/dev_data/dev_data.csv')

    from gzjwst.survey_products import download_mast_products_simplified
    download_mast_products_simplified(product_uris=df['dataURI'], save_dir='data/dev_data', max_files=None, cache=True)
