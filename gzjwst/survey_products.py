"""
Just a refactor of the below 'large download' script from the STScI GitHub repo.
https://github.com/spacetelescope/notebooks/blob/master/notebooks/MAST/Astroquery/large_downloads/companion_script.py
"""

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



def download_mast_products(products, save_dir, chunk_size=5, curl_flag=False):
    """
    Function that downloads products from MAST from a given table of observations

    Args:
        products (astropy table): product table as returned by get_observations_by_*
        save_dir (str): directory to save files in
        chunk_size (int, optional): Number of images to download per chunk if downloading large numbers of files. Defaults to 5 as suggested by the MAST team.
        curl_flag (bool, optional): If True, generates a curl download script instead of downloading data. Defaults to False.
    """

    if len(products) <= chunk_size:
        manifest = Observations.download_products(products, productType=['SCIENCE'], download_dir=save_dir, flat=True, curl_flag=curl_flag)


    else:
        print('Downloading '+str(len(products)) + ' files ('+str(round(products['size'].sum()/1e9,2))+' GB)')
        print('Splitting into chunks of '+str(chunk_size)+' files\n')
    
        chunks = [products[i:i+chunk_size] for i in range(0, len(products), chunk_size)]
        for i, chunk in enumerate(chunks):
            print('\nChunk '+str(i+1)+' of '+str(len(chunks)))
            print('Downloading '+str(len(chunk)) + ' files ('+str(round(chunk['size'].sum()/1e9,2))+' GB)\n')
            manifest = Observations.download_products(chunk, productType=['SCIENCE'], download_dir=save_dir, flat=True, curl_flag=curl_flag)

       

    #    TODO use requests or similar to download a few
    # or use curl_flag=True for script version, TBD
    # output should be level 3 products download to appropriate directory