"""
Just a refactor of the below 'large download' script from the STScI GitHub repo.
https://github.com/spacetelescope/notebooks/blob/master/notebooks/MAST/Astroquery/large_downloads/companion_script.py
"""
import time
import tqdm
import os
import pandas as pd
import logging
import asyncio

from astroquery.mast import Observations

def get_observations(**user_query_kwargs):
    """
    Get a table of observations for a given proposal ID or provenance name.
    query_kwargs could be: proposal_id=foo, provenance_name=bar, instrument_name=baz, etc.

    See https://astroquery.readthedocs.io/en/latest/api/astroquery.mast.ObservationsClass.html#astroquery.mast.ObservationsClass.query_criteria
    The Column Name is the keyword, with the argument being one or more acceptable values for that parameter, 
    except for fields with a float datatype where the argument should be in the form [minVal, maxVal]. 
    For non-float type criteria wildcards maybe used (both * and % are considered wildcards), 
    however only one wildcarded value can be processed per criterion
    """

    default_query_kwargs = {'instrument_name':'NIRCAM/IMAGE', 'dataRights': 'PUBLIC', 'intentType': 'SCIENCE'}
    query_kwargs = {**default_query_kwargs, **user_query_kwargs}
    # Get the list of observations for this proposal
    obs_table = Observations.query_criteria(**query_kwargs)
    logging.info(f'Found {len(obs_table)} observations, checking for products')

    # temp: debug
    # obs_table.to_pandas().to_csv('data/obs_table.csv', index=False)
    # exit()

    # temp: test the chunk size
    # start_time = time.time()
    # obs_table = obs_table[:30]
    # products = Observations.get_product_list(obs_table)
    # end_time = time.time()
    # logging.info(f'Got {len(products)} products in {end_time-start_time:.2f} seconds')
    # exit()

    products = get_product_list_batched(obs_table)
    return products


def get_product_list_batched(obs_table, chunk_size=30):
    # https://spacetelescope.github.io/mast_notebooks/notebooks/multi_mission/large_downloads/large_downloads.html#retreive-associated-products     
    
    # Split the observations into chunks and aysnc ask for associated products for each chunk
    chunks = [obs_table[i:i+chunk_size] for i in range(0, len(obs_table), chunk_size)]
    df = pd.concat([Observations.get_product_list(chunk).to_pandas() for chunk in tqdm.tqdm(chunks, unit='observation id chunk')])

    # Keep only the unique files
    df = df.drop_duplicates(subset='productFilename')

    logging.info(f'Found {len(df)} unique products')
    # will apply per-survey filters later, in survey script
    return df.reset_index(drop=True)


def download_mast_products_simplified(product_uris: list[str], save_dir: str, max_files=None, cache=True):
    """
    Dowload a list of products from MAST, using the URI returned by Observations.get_product_list.

    Args:
        product_uris (list[str]): list of uris to download
        save_dir (str): directory to save the files. Files named by URI within that directory
        max_files (int, optional): Download at most this many files, for debugging. Defaults to None.
        cache (bool, optional): Use MAST's size check to avoid redownloading existing files of the correct size. Defaults to True.
    """
    if (max_files is not None) and (len(product_uris) > max_files):
        logging.warning(f'Downloading only the first {max_files} of {len(product_uris)}')
        product_uris = product_uris[:max_files]
        
    for product_n, product_uri in enumerate(product_uris):
        logging.info(f'Downloading file {product_n+1} of {len(product_uris)}')
        filename = os.path.basename(product_uri)
        local_path = os.path.join(save_dir, filename)
        Observations.download_file(uri=product_uri, local_path=local_path, cache=cache)
