"""
Each cutout function should expect a catalog row
i.e.: a path to a (in general, multiple) fits files, plus coordinate and sizing metadata
Should also expect parameters controlling the cutout colouring, dynamic range, etc
"""


def make_single_band_cutout(galaxy: pd.Series, todo_more_params):
    pass  # should return nd array ready to be saved as jpg/png


# TODO etc for other band combinations, copy from Hayley's notebook


def select_2d_cutout():
    pass # TODO David or Hayley's slice code to select 2D flux values from a mosaic

