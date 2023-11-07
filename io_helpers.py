import os
from PIL import Image
from io import BytesIO
import requests
from importlib import reload
from glob import glob
from astropy.io import fits
from pyavm import AVM
from pyavm.exceptions import NoXMPPacketFound
import wcs_helpers as wh
reload(wh) 
import logger as logger
reload(logger)
FITS_EXTENSIONS = [".fits", ".fit", ".fts", ".fz", ".fits.fz"]
IMAGE_EXTENSIONS = [".jpg", ".jpeg", ".png"]



def _convertToSimpleJPEG(filename):
    """
    Open a JPEG file and save it again as a simple JPEG file.
    *** Warning: this will remove any AVM data from the file ***
    """
    im = Image.open(filename)
    im.save(filename)



def get_PIL_image(im):
    if isinstance(im, str):
        if im.startswith("http"):
            return get_PIL_image_from_url(im)
        # with Image.open(im) as pil_im:
        #     im = pil_im
        return Image.open(im)
    elif isinstance(im, Image):
        return im

def url_as_file(url):
    """
    Convert a URL to a filename
    """
    return BytesIO(requests.get(url).content)

def get_PIL_image_from_url(url):
    """
    Get a PIL image from a URL
    """
    return Image.open(url_as_file(url))



def get_image_size(image_path):
    """
    Get the image size from a PIL image
    
    returns (height, width)
    """
    # if it is a web file, return the size from the web
    if isinstance(image_path, str):
        if image_path.startswith("http"):
            im = get_PIL_image_from_url(image_path)
            return im.size
        with Image.open(image_path) as im:
            width, height = im.size
        return width, height
    
    elif isinstance(image_path, Image):
        return im.size
    elif isinstance(image_path, tuple):
        # assume it is a tuple of (width, height)
        return image_path

def get_suffix(suffix):
    # add "_" to the suffix if it doesn't start with "_" and isn't empty
    return ("_" if suffix and suffix.find("_") != 0 else "") + suffix

def open_header(filename):
    try:
        return fits.Header.fromfile(filename)
    except:
        return fits.Header.fromfile(filename, endcard=False, padding=False)
    

def get_image_header(image_file, wcs_file = None):
    """
    Get the image header from a FITS file.
    """
    if wcs_file is not None:
        return open_header(wcs_file)
    else:
        # check if the file is a FITS file
        if os.path.splitext(image_file)[1] in FITS_EXTENSIONS:
            logger.log("Reading FITS header from input FITS file", level="INFO")
            return open_header(image_file), None
        elif os.path.splitext(image_file)[1] in IMAGE_EXTENSIONS:
            try:
                logger.log("Reading AVM from image file", level="INFO")
                avm = AVM.from_image(image_file)
                return avm.to_wcs().to_header(), avm
            except:
                logger.log("Could not read AVM from image file", level="ERROR")
                logger.log("Returning blank header")
                return wh.blank_header(), None
                
def get_avm(filename):
    try:
        return AVM.from_image(filename)
    except NoXMPPacketFound:
         return None

def get_clean_header(wcsfile, remove_comments=True, remove_sip = True):
    if isinstance(wcsfile, str):
        logger.log(f"Reading header from file {wcsfile}", level='DEBUG')
    else:
        logger.log(f"io_helpers.get_clean_header: wcsfile is type{type(wcsfile)}", level='ERROR')
    header = open_header(wcsfile)
    if remove_comments:
        logger.log("Cleaning header", level='DEBUG')
        header = wh.clean_header(header)
    if remove_sip:
        logger.log("Removing SIP", level='DEBUG')
        header = wh.remove_sip(header)
    return header

def is_fits(filename):
    return os.path.splitext(filename)[1] in FITS_EXTENSIONS

def is_image(filename):
    return os.path.splitext(filename)[1] in IMAGE_EXTENSIONS
    
def guess_wcs_filename(image_path):
    """
    Guess the WCS filename from the image filename
    This should be the same name as the image but
    with .wcs or .wcs.fits or .fits extension
    it may also be a file called wcs.fits in the same directory
    """
    image_path = os.path.abspath(image_path)
    image_dir = os.path.dirname(image_path)
    image_name = os.path.basename(image_path)
    image_name = image_name.split(".")[0]
    wcs_files = glob(os.path.join(image_dir, f"{image_name}.wcs*"))
    if len(wcs_files) > 0:
        return wcs_files[0]
    else:
        wcs_files = glob(os.path.join(image_dir, "wcs.fits"))
        if len(wcs_files) > 0:
            return wcs_files[0]
        else:
            return None

def get_current_or_on_path(filename, path):
    """
    Get the current directory or on the path
    """
    if os.path.exists(filename):
        return filename
    else:
        if os.path.exists(os.path.join(path, filename)):
                return os.path.join(path, filename)
    return None
    

def github_raw_path(github_id, repository, branch, path):
    path = os.path.join(f"{github_id}/{repository}/{branch}", path)
    return "https://raw.githubusercontent.com/" + path


def convert_fits_to_jpg(filename):
    """
    reads in a fits file with 3 dimensions (color, x, y)
    and saves it as a jpg file
    """
    im = fits.open(filename)
    im = im[0].data
    # reshape the image to be (x, y, color)
    im = im.transpose(1, 2, 0)
    im = Image.fromarray(im)
    im.save(filename.replace(".fits", ".jpg"))
    