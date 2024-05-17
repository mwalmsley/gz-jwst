if __name__ == "__main__":

    from astroquery.mast import Observations


    meta_table = Observations.get_metadata("observations")
    print(meta_table.to_pandas())