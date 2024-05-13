

def get_filenames_from_hlsp(hlsp_name, save_dir):
    Observations.query_criteria(provenance_name=hlsp_name)
    # and save the filenames to relevant directory