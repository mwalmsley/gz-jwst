# gz-jwst

## MASTer plan

JWST morphology science has mostly been *a* team looking at *their* survey.
We would like do science across many JWST surveys.
This will help with scale (which we're good at) and will lead to results that unify the otherwise-incomparable measurements of e.g. bar fraction.

We will use the calibrated data products available via MAST. This loses some possible improvements from individual teams but ensures processing consistency between surveys.

We will load these data products via DataLabs which avoids needing to download them.

We will make our own catalogs with a simple SEP script aiming to do a "good enough" job. Our target sources are big and easy to find. Using our own catalog gives us control of the selection function.

We will then make cutouts based on our catalogs, though we're still a bit unsure of the best way to colour these and present them to volunteers.

## Next Steps

We already have the code to do most if not all of this. We're putting it on GitHub and sticking it together.

1. Hayley's astroquery.mast script for identifying observation_id in each survey
2. David's lookup script to find mosiacs for the observation_id within the DataLabs volumes
3. Kameswara's SEP script to make catalogs
4. David's cutout script to extract single-band cutouts from the mosaics
5. Hayley's colouring to make images (probably with tweaks)
