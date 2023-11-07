from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

# we want to create a class that takes in a fits header
# this class should have a property that returns a WCS object
# we should be able to rotate the header by an angle
# we should be able to shift the header by a given amount
# we should be able to scale the header by a given amount
# the header should be normalized to a valid FITS header
# using CDELT and CD matrix (not the PC matrix which astropy uses)
