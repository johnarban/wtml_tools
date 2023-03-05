import os
from PIL import Image
from importlib import reload
from glob import glob
from astropy.io import fits
from pyavm import AVM
from pyavm.exceptions import NoXMPPacketFound
import wcs_helpers as wh
reload(wh) 
from logger import log, set_debug_level
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
        with Image.open(im) as pil_im:
            im = pil_im
        return im
    elif isinstance(im, Image):
        return im

def get_image_size(image_path):
    """
    Get the image size from a PIL image
    
    returns (height, width)
    """
    if isinstance(image_path, str):
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

def get_image_header(image_file, wcs_file = None):
    """
    Get the image header from a FITS file.
    """
    if wcs_file is not None:
        return fits.Header.fromfile(wcs_file)
    else:
        # check if the file is a FITS file
        if os.path.splitext(image_file)[1] in FITS_EXTENSIONS:
            log("Reading FITS header from image file", level=1)
            return fits.Header.fromfile(image_file), None
        elif os.path.splitext(image_file)[1] in IMAGE_EXTENSIONS:
            try:
                log("Reading AVM from image file", level=1)
                avm = AVM.from_image(image_file)
                return avm.to_wcs().to_header(), avm
            except:
                log("Could not read AVM from image file", level=2)
                log("Returning blank header")
                return wh.blank_header(), None
                
def get_avm(filename):
    try:
        return AVM.from_image(filename)
    except NoXMPPacketFound:
         return None

def get_clean_header(wcsfile, remove_comments=True, remove_sip = True):
    if isinstance(wcsfile, str):
        log(f"Reading header from file {wcsfile}", level='INFO')
    else:
        print(wcsfile)
    header = fits.Header.fromfile(wcsfile)
    if remove_comments:
        log("Cleaning header", level='DEBUG')
        header = wh.clean_header(header)
    if remove_sip:
        log("Removing SIP", level='DEBUG')
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

