from  aplpy import FITSFigure

from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u 
from astropy.wcs import WCS

from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
import wcs_helpers as wh
import numpy as np


DEBUG_LEVELS = {"DEBUG": 0, "INFO": 1}
DEBUG_LEVEL = 1
FITS_EXTENSIONS = [".fits", ".fit", ".fts", ".fz", ".fits.fz"]
IMAGE_EXTENSIONS = [".jpg", ".jpeg", ".png"]


def set_debug_level(level):
    global DEBUG_LEVEL
    if isinstance(level, str):
        DEBUG_LEVEL = DEBUG_LEVELS[level]
    else:
        DEBUG_LEVEL = level


def log(arg: str, level=0):
    # check if level is a string
    if isinstance(level, str):
        level = DEBUG_LEVELS[level]

    if level >= DEBUG_LEVEL:
        print(arg)



def preview_image(url,coords=None):
    g =FITSFigure(url)
    g.add_grid()
    g.set_theme('pretty')
    
    if coords is None:
        xy = g._data.shape # get coords at center
        # get coords at center
        ra,dec = g.pixel2world(min(*xy),min(*xy))
        coord_center = SkyCoord(ra=ra,dec=dec,unit=(u.deg, u.deg))
        catalog = Vizier.query_region(coord_center,radius=2*u.deg,catalog='III/135A/catalog') # HD catalog
        coords = catalog[0]['_RA.icrs','_DE.icrs']
        coords.rename_columns(['_RA.icrs','_DE.icrs'],['ra','dec'])
        coords = SkyCoord(coords['ra'],coords['dec'],unit=(u.hourangle, u.deg))

    g.show_markers(coords.ra,coords.dec,marker='x',facecolor='red',edgecolor='red',s=15)
    g.axis_labels.show()
    # g.tick_labels.set_xformat('hh:mm:ss.ss')
    g.tick_labels.show()
    g.tick_labels.set_style('colons')
    g.show_rgb()
    return g
    

def rgb_reproject(image, header):
    from_wcs = WCS(header)
    r, g, b = image.split()
    r = np.array(r)[:, ::-1]
    g = np.array(g)[:, ::-1]
    b = np.array(b)[:, ::-1]
    to_wcs, shape_out = find_optimal_celestial_wcs([(r.shape, from_wcs)])
    r, _ = reproject_interp((r, from_wcs), to_wcs, shape_out=shape_out)
    g, _ = reproject_interp((g, from_wcs), to_wcs, shape_out=shape_out)
    b, _ = reproject_interp((b, from_wcs), to_wcs, shape_out=shape_out)
    out = np.nan_to_num([r, g, b])
    return out, to_wcs




def do_parity_inversion(header, im, force=False):
    # flip image top to bottom
    # from https://github.com/WorldWideTelescope/wwt-aligner/blob/master/backend/wwt_aligner/driver.py#L449
    parity = wh.get_parity(header=header)
    if (parity < 0) and not force:
        log(
            "Not inverting parity: WWT needs JPEGs to have negative parity",
            level="DEBUG",
        )
        # WWT requires images with negative parity
        return header
    log("Inverting parity", level="INFO")
    # wcs=WCS(header)
    # hdwork=wcs.to_header(relax=False)
    hdwork = wh.flip_parity(header, im.height)
    if "NAXIS1" not in header:
        hdwork = wh.add_NAXES(hdwork, *im.size[::-1])
    return hdwork


def add_scale_rot(header, from_header=None):
    header = header.copy()
    if from_header is None:
        scale, rotation, parity = wh.get_scale_rot(header)
    else:
        scale, rotation, parity = wh.get_scale_rot(from_header)
    # if parity > 0:
    #     rotation=180 - rotation
    header["CDELT1"] = scale[0]
    header["CDELT2"] = scale[1]
    header["CROTA2"] = rotation
    # remove CD matrix
    # header['CRPIX2']=header['NAXIS2'] + 1 - header['CRPIX2']
    return header