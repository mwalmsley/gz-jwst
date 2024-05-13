"""
Just a refactor of the below 'large download' script from the STScI GitHub repo.
https://github.com/spacetelescope/notebooks/blob/master/notebooks/MAST/Astroquery/large_downloads/companion_script.py
"""


def get_filenames_for_proposal(proposal_id: str) -> List[str]:
    """
    Get the filenames for a given proposal ID.
    """
    # Get the list of observations for this proposal
    obs_table = Observations.query_criteria(proposal_id=proposal_id)
    obs_ids = obs_table['obs_id']

    # Get the list of products for each observation
    products = Observations.get_product_list(obs_ids)
    # Get the list of filenames for each product
    filenames = Observations.get_filenames(products)
    return filenames

def download_filenames(filenames, save_dir, max_files=3):
       manifest = Observations.download_products(
               files
               , curl_flag=False
               , productType=['SCIENCE']
               )
    #    TODO use requests or similar to download a few
    # or use curl_flag=True for script version, TBD
    # output should be level 3 products download to appropriate directory