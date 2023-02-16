import os

# from lxml.etree import ElementTree as ET
# import lxml.etree as etree
import xml.etree as etree
import xml.etree.ElementTree as ET
from glob import glob
from importlib import reload
from io import BytesIO

import numpy as np
from astropy.io.fits import Header
from astropy.wcs import WCS
from astropy.wcs.utils import _is_cd_orthogonal, is_proj_plane_distorted
from PIL import Image
from pyavm import AVM, NoAVMPresent
from pyavm.exceptions import NoXMPPacketFound
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs

import wcs_helpers as wh

reload(wh)
from xml_indenter import smart_indent_xml, split_xml_attributes

DEBUG_LEVELS = {"DEBUG": 0, "INFO": 1}
DEBUG_LEVEL = 1


def set_debug_level(level):
    global DEBUG_LEVEL
    if isinstance(level, str):
        DEBUG_LEVEL = DEBUG_LEVELS[level]
    else:
        DEBUG_LEVEL = level


def log(arg, level=None):
    # check if level is a string
    if isinstance(level, str):
        level = DEBUG_LEVELS[level]

    if level >= DEBUG_LEVEL:
        print(arg)


def _convertToSimpleJPEG(filename):
    """
    Open a JPEG file and save it again as a simple JPEG file.
    *** Warning: this will remove any AVM data from the file ***
    """
    im = Image.open(filename)
    im.save(filename)


def get_header(wcsfile, clean=False, clean_sip=False, clean_all=False):
    header = Header.fromfile(wcsfile)
    if clean or clean_all:
        log("Cleaning header", level='DEBUG')
        header = wh.clean_header(header)
    if clean_sip or clean_all:
        log("Removing SIP", level='DEBUG')
        header = wh.remove_sip(header)
    return header


def get_suffix(suffix):
    return ("_" if suffix.find("_") != 0 else "") + suffix

def get_PIL_image(im):
    if isinstance(im, str):
        return Image.open(im)
    elif isinstance(im, Image):
        return im

def output_avm(im, header, name="image", suffix="", path_out=".", ext="jpg", remove_full_fits_header=False):
    """
    output_avm takes an image and a header and embeds the AVM in the image
    write's the image out to path_out/name_suffix_tagged.ext

    Parameters
    ----------
    im : str or PIL.Image
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
        extension of output image, by default "jpg"
    remove_full_fits_header : bool, optional
        remove full FITS header from AVM, by default False

    Returns
    -------
    str, pyavm.AVM
        path to image tagged with AVM, AVM object
        defaul output is ./image_tagged.jpg
    """
    
    # ensure we have a PIL image
    im = get_PIL_image(im)
    
    log(f"\n ****** Embedding AVM in {ext} file ****** \n", level='INFO')
    name = name + get_suffix(suffix)

    save_image_path = os.path.join(path_out, name + "_tagged." + ext)
    log(f"AVM tagged image saved to {save_image_path}", level='INFO')

    if isinstance(header, WCS):
        header = header.to_header(relax=True)
        header = wh.add_NAXES(header, im)

    hd = wh.make_avm_header(header)  # puts scale, rot in header, removing cd matrix
    avm = AVM.from_header(hd)
    if remove_full_fits_header: avm.Spatial.FITSheader = ""
    log(avm.Spatial, level='DEBUG')
    
    
    temp = "temp_can_delete." + ext 
    im.save(temp)
    avm.embed(temp, save_image_path)
    os.remove(temp)

    return save_image_path, avm


def add_avm_tags(image_path, wcsfile="wcs.fits", name=None, path_out=".", suffix=""):
    im = Image.open(image_path)
    ext = image_path.split(".")[-1]
    suffix = get_suffix(suffix)

    # TODO: preferred behavior
    # if image is FITS file, use astropy to get header,
    # if has AVM and no WCS then use AVM
    # always separate WCS file if it is specified or exists
    if wcsfile is None:
        wcsfile = image_path.replace("." + ext, ".wcs")
    elif wcsfile.find('*') == 0:
        wcsfile = glob(wcsfile)[0]

    if name is None:
        name = image_path.replace(f".{ext}", "")

    header = get_header(wcsfile, clean_all=True)
    header = wh.add_NAXES(header, im, add_naxisi=True)

    out_tagged, avm = output_avm(im, header, name=name, suffix=suffix, path_out=path_out, ext=ext)
    return header, WCS(header), avm, out_tagged, im


def wtml_header(header):
    header = header.copy()
    has_cd = "CD1_1" in header
    if has_cd:
        log("\tHeader has CD matrix", level='DEBUG')
    has_crota = "CROTA2" in header
    if has_crota:
        log("\tHeader has CROTA", level='DEBUG')
    has_cdelt = "CDELT1" in header
    if has_cdelt:
        log("\tHeader has CDELT", level='DEBUG')
    has_pc = "PC1_1" in header
    if has_pc:
        log("\tHeader has PC matrix", level='DEBUG')

    wcs = WCS(header)
    crpix = [(header["NAXIS1"] + 1) / 2, (header["NAXIS2"] + 1) / 2]
    crval = wcs.wcs_pix2world(*crpix, 1)
    offset = crpix
    cd = wh.get_cd(wcs)
    scales, rot, parity = wh.get_scale_rot(header)
    scale = np.sqrt(np.abs(scales[0] * scales[1]))

    return rot, crval, offset, scale, cd, parity


def create_wtml(header, im, name="test", url=None, credits="credits", credits_url="credits_url", thumbnail_url=None, out=None):
    # print out header sayin "Creating WTML file"
    log("\n ****** Creating WTML file ****** \n", level='INFO')
    if out is None:
        out = name + ".wtml"
    log("WTML file:", out)

    if url is None:
        raise ValueError("I need a url")

    rot, crval, offset, scale, cd, parity = wtml_header(header)

    imageset = create_imageset_dict(name, rot, crval, offset, scale, im.height, im.width, url, parity, cd)

    folder = set_folder(name)
    place = set_place(name, crval[0] / 15, crval[1])

    # log('Creating WTML file')
    el_imageset = ET.Element("ImageSet", attrib=imageset)

    el_credits = ET.Element("Credits")
    el_credits.text = credits

    el_credits_url = ET.Element("CreditsUrl")
    el_credits_url.text = credits_url

    el_description = ET.Element("Description")
    el_description.text = "description"

    el_thumbnail_url = ET.Element("ThumbnailUrl")
    el_thumbnail_url.text = thumbnail_url or "https://nova.astrometry.net/image/16942765"

    el_imageset.append(el_description)
    el_imageset.append(el_credits)
    el_imageset.append(el_credits_url)
    el_imageset.append(el_thumbnail_url)

    el_forground = ET.Element("ForegroundImageSet")
    el_forground.append(el_imageset)

    el_place = ET.Element("Place", attrib=place)
    el_place.append(el_forground)

    el_folder = ET.Element("Folder", attrib=folder)
    el_folder.append(el_place)

    smart_indent_xml(el_folder)
    tree = ET.ElementTree(el_folder)
    tree.write(out)
    split_xml_attributes(out, field="ImageSet")
    return tree, imageset


def create_imageset_dict(name, rot, crval, offset, scale, height, width, url, parity, cd):
    log("\n WTML Place & Imageset values: \n", level=3)
    # log('  name', name)
    # log('  parity', parity)
    bottoms_up = parity > 0
    # log('  bottoms_up', bottoms_up)
    # log('  rot', rot)
    # log('  crval', crval)
    # log('  offset', offset)
    # log('  scale', scale)
    # log('  height', height)
    # log('  width', width)
    full_string = f"Name: {name} Rotation: {rot:0.2f} Bottoms-up: {bottoms_up} Center_X: {crval[0]:0.2f} Center_Y: {crval[1]:0.2f} Offset_X: {offset[0]} Offset_Y: {offset[0]} Scale: {scale:0.3g} Height: {height} Width: {width} URL: {url}"
    log(full_string, level=2)
    if "github" in url:
        # repository is the github repo in the url
        repository = "/".join(url.split("/")[3:5])
        log(f"  repository @{repository}")
        directory = "/".join(url.split("/")[4:-1])
        log(f"  directory: {directory}/")
    file = url.split("/")[-1]
    log(f"\timage: {file}")

    log(f"\turl: {url}\n")

    imageset = {
        "DemUrl": "",
        "MSRCommunityId": "0",
        "MSRComponentId": "0",
        "Permission": "0",
        "Generic": "False",
        "DataSetType": "Sky",
        "BandPass": "Visible",
        "Url": url,
        "TileLevels": "0",
        "WidthFactor": "2",
        "Sparse": "True",
        "Rotation": str(rot),
        "QuadTreeMap": "",
        "Projection": "SkyImage",
        "Name": name,
        "FileType": ".jpg",
        "CenterY": str(crval[1]),  # ra in decimal degrees
        "CenterX": str(crval[0]),  # dec in decimal degrees
        "BottomsUp": str(bottoms_up),  # True if it should be inverted
        "StockSet": "False",
        "ElevationModel": "False",
        "OffsetX": str((offset[0])),  # ra offset in pixels
        "OffsetY": str((offset[1])),  # dec offset in pixels
        "BaseTileLevel": "0",
        "BaseDegreesPerTile": "",
        "ReferenceFrame": "",
        "MeanRadius": "0",
    }

    if imageset["Projection"] == "SkyImage":
        imageset["BaseDegreesPerTile"] = str(scale)  # str(np.sqrt(cd[0, 1]**2 + cd[1, 1]**2)) # scal_y
    else:
        imageset["BaseDegreesPerTile"] = str(max(width, height) * scale)  # str(np.sqrt(cd[0, 1]**2 + cd[1, 1]**2)) #str(height_arcmin * max(height, width) / (60 * height))

    base_url = "http://www.worldwidetelescope.org/wwtweb/ShowImage.aspx?"
    options = (bottoms_up, 3600 * scale, name, url, crval[0], crval[1], offset[0], offset[1], rot + 180, "https://nova.astrometry.net/image/16942765")
    wwt_url = "reverseparity=%s&scale=%.6f&name=%s&imageurl=%s&ra=%.6f&dec=%.6f&x=%.1f&y=%.1f&rotation=%.2f&thumb=%s" % options
    log("WorldWide Telescope URL:")
    log(base_url + wwt_url, level=2)
    log("\n", level=3)

    return imageset


def set_folder(name):
    folder = {
        "MSRCommunityId": "0",
        "MSRComponentId": "0",
        "Permission": "0",
        "Name": name,
        "Group": "Explorer",
        "Searchable": "False",
    }
    return folder


def set_place(name, ra, dec, rotation=0, zoom=10):
    place = {
        "Name": name,
        "DataSetType": "Sky",
        "RA": str(ra),  # ra in hours
        "Dec": str(dec),  # dec in decimal degrees
        "Rotation": str(rotation),
        "ZoomLevel": str(zoom),
    }
    return place


def parse_xml(xml):
    try:
        return ET.parse(xml)
    except:
        parser = etree.XMLParser(remove_blank_text=True)
        return etree.parse(xml, parser)
    return None


def set_xml_element(el, attribute=None, text=None):
    if attribute is not None:
        try:
            el.attrib = attribute
        except:
            for key, value in attribute.items():
                el.set(key, value)

    if text is not None:
        el.text = text


def rgb_reproject(image, header):
    from_wcs = WCS(header)
    r, g, b = image.split()
    r = np.array(r)[:, ::-1]
    g = np.array(g)[:, ::-1]
    b = np.array(b)[:, ::-1]
    to_wcs, shape_out = find_optimal_celestial_wcs([(r.shape, from_wcs)])
    r, footprint = reproject_interp((r, from_wcs), to_wcs, shape_out=shape_out)
    g, footprint = reproject_interp((g, from_wcs), to_wcs, shape_out=shape_out)
    b, footprint = reproject_interp((b, from_wcs), to_wcs, shape_out=shape_out)
    out = np.nan_to_num([r, g, b])
    return out, to_wcs


def avm_to_wtml(avm):
    avm_header = avm.to_wcs().to_header()
    avm_header["NAXIS"] = 2
    avm_header["NAXIS1"] = avm.Spatial.ReferenceDimension[0]
    avm_header["NAXIS2"] = avm.Spatial.ReferenceDimension[1]

    avm_header = wh.flip_parity(avm_header)

    return wtml_header(avm_header)


def add_scale_rot(header, from_header=None, wwt=False):
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


def do_parity_inversion(header, im, force=False):
    # flip image top to bottom
    # from https://github.com/WorldWideTelescope/wwt-aligner/blob/master/backend/wwt_aligner/driver.py#L449
    parity = wh.get_parity(header=header)
    if (parity < 0) and not force:
        log("Not inverting parity: WWT needs JPEGs to have negative parity", level='DEBUG')
        # WWT requires images with negative parity
        return header
    log("Inverting parity", level='INFO')
    # wcs=WCS(header)
    # hdwork=wcs.to_header(relax=False)
    hdwork = wh.flip_parity(header, im.height)
    if "NAXIS1" not in header:
        hdwork = wh.add_NAXES(hdwork, im)
    return hdwork
