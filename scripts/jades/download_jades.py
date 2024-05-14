from gzjwst.survey_products import get_observations, download_mast_products

PROVENANCE_NAME = 'jades'
DATA_DIRECTORY = 'data/jades/hlsp'
# NIRCAM_FILTERS = ['F090W', 'F115W', 'F150W', 'F182M', 'F200W', 'F210M', 'F277W', 'F335M', 'F356W', 'F410M', 'F444W']

# get all data products associated with CEERS
data_products = get_observations(provenance_name=PROVENANCE_NAME, query_kwargs={'instrument_name':'NIRCAM/IMAGE'})


# now split into goods-s and goods-n data
goods_s = [i for i in range(len(data_products)) if 'goods-s' in data_products['obs_id'][i]]
goods_n = [i for i in range(len(data_products)) if 'goods-n' in data_products['obs_id'][i]]

data_products_goods_s = data_products[goods_s]
data_products_goods_n = data_products[goods_n]


# test call - downloads just a script
download_mast_products(data_products_goods_s, save_dir=DATA_DIRECTORY+'goods_s/', curl_flag=True)
download_mast_products(data_products_goods_n, save_dir=DATA_DIRECTORY+'goods_n/', curl_flag=True)

# real script
# download_mast_products(data_products_goods_s, save_dir=DATA_DIRECTORY+'goods_s/')
# download_mast_products(data_products_goods_n, save_dir=DATA_DIRECTORY+'goods_n/')