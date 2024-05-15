"""
Just a refactor of the below 'large download' script from the STScI GitHub repo.
https://github.com/spacetelescope/notebooks/blob/master/notebooks/MAST/Astroquery/large_downloads/companion_script.py
"""
import os
import logging

from astroquery.mast import Observations

def get_observations_by_proposal_id(proposal_id: str):
    """
    Get a table of observations for a given proposal ID.
    """
    # Get the list of observations for this proposal
    obs_table = Observations.query_criteria(proposal_id=proposal_id)

    # Get the list of products for each observation
    products = Observations.get_product_list(obs_table)

    return products

def get_observations_by_provenance_name(provenance_name: str):
    """
    Get a table of observations for a given provenance name.
    """
    # Get the list of observations for this proposal
    obs_table = Observations.query_criteria(provenance_name=provenance_name)

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



# def download_mast_products(products, save_dir, max_files=3, chunk_size=None, curl_flag=False):
#     """
#     Function that downloads products from MAST from a given table of observations

#     Args:
#         products (astropy table): product table as returned by get_observations_by_*
#         save_dir (str): directory to save files in
#         max_files (int, optional): The maximum number of files to download. Defaults to 3.
#         chunk_size (int, optional): Number of images to download per chunk if downloading large numbers of files. Defaults to None. Max is 5.
#         curl_flag (bool, optional): If True, generates a curl download script instead of downloading data. Defaults to False.
#     """

#     if len(products) > max_files and not chunk_size:
#         products = products[:max_files]
#         print('More filenames provided than the maximum number of files to download')
#         print('(Set chunk_size if you wish to download files in chunks)')
#         print('Downloading the first '+str(max_files)+' files('+str(round(products['size'].sum()/1e9,2))+' GB)\n')

#         manifest = Observations.download_products(products, productType=['SCIENCE'], download_dir=save_dir, flat=True, curl_flag=curl_flag)

#     elif len(products) <= max_files and not chunk_size:
#         print('Downloading '+str(len(products)) + ' files ('+str(round(products['size'].sum()/1e9,2))+' GB)\n')
#         manifest = Observations.download_products(products, productType=['SCIENCE'], download_dir=save_dir, flat=True, curl_flag=curl_flag)


#     elif len(products) <= max_files and chunk_size:
#         if chunk_size > 5:
#             chunk_size = 5
#             print('Chunk size cannot be bigger than 5. Defaulting to 5.\n')

#         print('Downloading '+str(len(products)) + ' files ('+str(round(products['size'].sum()/1e9,2))+' GB)')
#         print('Splitting into chunks of '+str(chunk_size)+' files\n')
    
#         chunks = [products[i:i+chunk_size] for i in range(0, len(products), chunk_size)]
#         for i, chunk in enumerate(chunks):
#             print('\nChunk '+str(i+1)+' of '+str(len(chunks)))
#             print('Downloading '+str(len(chunk)) + ' files ('+str(round(chunk['size'].sum()/1e9,2))+' GB)\n')
#             manifest = Observations.download_products(chunk, productType=['SCIENCE'], download_dir=save_dir, flat=True, curl_flag=curl_flag)

       

    #    TODO use requests or similar to download a few
    # or use curl_flag=True for script version, TBD
    # output should be level 3 products download to appropriate directory

# TODO replace the curl scripts with our own wget/curl that uses the list of urls from file option, pointed at the MAST urls like these
"""
curl --globoff --location-trusted -f --progress-bar --create-dirs $CONT --output ./${DOWNLOAD_FOLDER}'/HLSP/hlsp_ceers_jwst_nircam_nircam2_f200w_v0.5_i2d/hlsp_ceers_jwst_nircam_nircam2_f200w_v0.5_i2d.fits.gz' 'https://mast.stsci.edu/api/v0.1/Download/file?bundle_name=MAST_2024-05-14T1721.sh&uri=mast:HLSP/ceers/nircam/hlsp_ceers_jwst_nircam_nircam2_f200w_v0.5_i2d.fits.gz'
"""