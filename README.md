# gz-jwst

## MASTer plan

JWST morphology science has mostly been *a* team looking at *their* survey.
We would like do science across many JWST surveys.
This will help with scale (which we're good at) and will lead to results that unify the otherwise-incomparable measurements of e.g. bar fraction.

We will use the calibrated data products available via MAST. This loses some possible improvements from individual teams but ensures processing consistency between surveys.

We will load these data products via DataLabs which avoids needing to download them.

We will make our own catalogs with a simple SEP script aiming to do a "good enough" job. Our target sources are big and easy to find. Using our own catalog gives us control of the selection function.

We will then make cutouts based on our catalogs, though we're still a bit unsure of the best way to colour these and present them to volunteers.

## What we're sticking together

We already have the code to do most if not all of this. We're putting it on GitHub and sticking it together.

1. Hayley's astroquery.mast script for identifying observation_id in each survey
2. David's lookup script to find mosiacs for the observation_id within the DataLabs volumes
3. Kameswara's SEP script to make catalogs
4. David's cutout script to extract single-band cutouts from the mosaics
5. Hayley's colouring to make images (probably with tweaks)

## Current Status

- Download FITS from MAST for all surveys &#10003;
- Find existing FITS on Datalabs &#10003; (but some are missing)
- Create catalog from FITS: WIP, works on one CEERS band but decisions needed to generalise
- Create cutouts from catalog/FITS: WIP, tested several options but decisions needed on final choice(s)

## Install

    conda create --name gz-jwst python>=3.9
    conda activate gz-jwst
    pip install -e .
    pip install -U --pre astroquery

Astroquery is separate because it has a [weird pip release style](https://astroquery.readthedocs.io/en/latest/).

## Code Structure

- Each survey has scripts under `/scripts` for downloading, creating a catalog, and creating cutouts.
- Those scripts take action by calling the same shared code (the `/gzjwst` package) with different parameters.
- We can experiment and pick good parameters using the playground notebooks in `notebooks`.

## Running

For each survey, run the `scripts/{survey}/download_{survey}.py` script. This first creates a list of products (FITS tiles) to download, and then downloads them from MAST.



## How the Pipeline Runs for All Surveys

### Look up all `observation_id` and get filenames

Find all `observation_id` for survey:

- For surveys without HLSP: Use astroquery.mast to get `observation_id` for all nircam observations that corresponds to a given proposal id. See script linked on notes. For each observation_id, get the filenames of all associated data products.
- For surveys with HLSP use astroquery.mast to pull the HLSP e.g. `Observations.query_criteria(provenance_name="ceers")`. Already returns productFilename.

For CEERS, we only have an HLSP version. For JADES, we have both LE3 and HLSP. For COSMOS, we have LE3. For PRIMER, we have nothing yet. For all surveys, we prefer LE3 for consistency!

For now: download a few and save locally. See `data` dir.


### Create cutouts from the catalog

Load a catalog table - either ASTRODEEP or COSMOS Webb's internal catalog. Slice out, from each corresponding mosaic, a 2D array of flux values around each source (with a variable field of view depending on the source size).

## Data

The most likely first surveys, selected simply by size, are:

- CEERS (including public data release v0.6, 200arcmin2)
- JADES (190arcmin2, possibly including GOODS-N as well as GOODS-N if it's on MAST yet)
- COSMOS-WEB (77arcmin2 now but 900arcmin2 from August)
- PRIMER TODO, and all other astrodeep surveys 
