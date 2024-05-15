"""
Just a refactor of the below 'large download' script from the STScI GitHub repo.
https://github.com/spacetelescope/notebooks/blob/master/notebooks/MAST/Astroquery/large_downloads/companion_script.py
"""
import os
import logging

from astroquery.mast import Observations

def get_observations(provenance_name=None, proposal_id=None, query_kwargs={}):
    """
    Get a table of observations for a given proposal ID or provenance name.
    """
    # Get the list of observations for this proposal
    if proposal_id and not provenance_name:
        obs_table = Observations.query_criteria(proposal_id=proposal_id, **query_kwargs)
    elif provenance_name and not proposal_id:
        obs_table = Observations.query_criteria(provenance_name=provenance_name, **query_kwargs)
    else:
        print('You must provide either a proposal ID or provenance name, but not both.')
        return

    # Get the list of products for each observation
    products = Observations.get_product_list(obs_table)

    return products


def download_mast_products_simplified(products, save_dir, max_files=None, cache=True):
    if (max_files is not None) and (len(products) > max_files):
        logging.warning(f'Downloading only the first {max_files} of {len(products)}')
        products = products[:max_files]
        
    for product_n, product_uri in enumerate(products['dataURI']):
        logging.info(f'Downloading file {product_n+1} of {len(products)}')
        filename = os.path.basename(product_uri)
        local_path = os.path.join(save_dir, filename)
        Observations.download_file(uri=product_uri, local_path=local_path, cache=cache)
