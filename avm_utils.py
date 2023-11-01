import os
from glob import glob
from importlib import reload
from PIL import Image
from pyavm import AVM
from astropy.io import fits
from astropy.io.fits import Header
from astropy.wcs import WCS
import wcs_helpers as wh
reload(wh)
import io_helpers as ih
reload(ih)
import helper_classes as hc

DEBUG_LEVELS = {"DEBUG": 0, "INFO": 1}
DEBUG_LEVEL = 1
FITS_EXTENSIONS = [".fits", ".fit", ".fts", ".fz", ".fits.fz"]
IMAGE_EXTENSIONS = [".jpg", ".jpeg", ".png"]

import logger as logger
reload(logger)

def make_avm_header(header):
    """
    Make an AVM header from a FITS header.
    Performs a parity flip if the parity is JPEG-like.
    """
    header = header.copy()
    header = wh.add_cd_matrix(header)

    parity = wh.get_parity(header=header)
    flip_parity = False
    if wh.is_JPEGLike(parity):
        flip_parity = True
        logger.log("This has JPEG-like parity", level=0)
        logger.log("flip parity for avm", level=1)
        header = wh.flip_parity(header)

    scale, rot, parity = wh.get_scale_rot(header)

    cdelt1 = scale[0]
    cdelt2 = scale[1]
    crota = rot

    if flip_parity:  # not wh.is_JPEGLike(parity):
        logger.log("Parity was flipped so, flipping cdelt2 (lat axis)", level=1)
        cdelt2 = -cdelt2
        cdelt1 = cdelt1
        crota = crota

    header["CDELT1"] = cdelt1
    header["CDELT2"] = cdelt2
    header["CROTA2"] = crota
    header = wh.remove_cd(header)
    return header



def write_avm(image, header, name="image", suffix="", path_out=".", ext="jpg", remove_full_fits_header=False):
    """
    output_avm takes an image and a header and embeds the AVM in the image
    write's the image out to path_out/name_suffix_tagged.ext

    Parameters
    ----------
    im : image path or PIL.Image
        image to embed AVM in
    header : astropy.io.fits.Header or astropy.wcs.WCS
        header to embed in AVM
    name : str, optional
        name of image, by default "image"
    suffix : str, optional
        suffix to add to image name, by default ""
    path_out : str, optional
        path to save image, by default "."
    ext : str, optional
        extension of output image, by default  same as input or "jpg" if PIL.Image
    remove_full_fits_header : bool, optional
        remove full FITS header from AVM, by default False

    Returns
    -------
    str, pyavm.AVM
        path to image tagged with AVM, AVM object
        defaul output is ./image_tagged.jpg
    """
    
    if isinstance(image, str):
        ext = image.split(".")[-1]
    
    if isinstance(image, str) & (name == "image"):
        name = os.path.basename(image).split(".")[0]
        if path_out == ".":
            path_out = os.path.dirname(image)
    
    # ensure we have a PIL image
    # image = ih.get_PIL_image(image)
    
    
    logger.log(f"\n ****** Embedding AVM in {ext} file ****** \n", level='INFO')
    name = name + ih.get_suffix(suffix)

    save_image_path = os.path.join(path_out, name + "_tagged." + ext)

    if header is None:
        if isinstance(image, str):
            ih.guess_wcs_filename(image_path = image)
            header = hc.ImageHeader(image).header


    # make_avm_header puts scale, rot in header, removes cd matrix
    # and applies parity flip if needed
    # header = make_avm_header(header)
    logger.log(header.__repr__(),level="DEBUG")
    avm = AVM.from_header(header)
    if remove_full_fits_header: avm.Spatial.FITSheader = ""
    avm.Spatial.Rotation = avm.Spatial.Rotation - 180
    logger.log(avm.Spatial, level='DEBUG')
    
    
    temp = "temp_can_delete." + ext 
    with Image.open(image) as im:
        im.save(temp)
    avm.embed(temp, save_image_path)
    # os.remove(temp)
    
    logger.log(f"AVM tagged image saved to {save_image_path}", level='INFO')

    return save_image_path, avm




def add_avm_tags(image_path, wcsfile=None, name=None, path_out=".", suffix=""):
    logger.log(f"\n ****** Processing {image_path} ****** \n", level='INFO')
    im = Image.open(image_path)
    ext = image_path.split(".")[-1]
    suffix = ih.get_suffix(suffix)

    if wcsfile is None:
        wcsfile = image_path.replace("." + ext, ".wcs")
    elif wcsfile.find('*') == 0:
        wcsfile = glob(wcsfile)[0]

    if name is None:
        name = image_path.replace(f".{ext}", "")
    
    image_header = hc.ImageHeader(image_path, wcsfile=wcsfile)
    header =image_header.header
    
    header = ih.get_clean_header(wcsfile)
    header = wh.add_NAXES(header, *im.size)

    out_tagged, avm = write_avm(im, header, name=name, suffix=suffix, path_out=path_out, ext=ext)
    return header, WCS(header), avm, out_tagged, im

