import os
from tokenize import String

# from lxml.etree import ElementTree as ET
# import lxml.etree as etree
import xml.etree as etree
import xml.etree.ElementTree as ET
from glob import glob
from importlib import reload
from io import BytesIO

import numpy as np
from astropy.io import fits
from astropy.io.fits import Header
from astropy.wcs import WCS
from PIL import Image

from pyavm import AVM

import wcs_helpers as wh
import io_helpers as ih
import avm_utils as au
import helper_classes as hc
import path_helpers as ph

reload(wh)
reload(ih)
reload(au)
reload(hc)
reload(ph)
import xml_indenter as xml
reload(xml)

from logger import log, set_debug_level
FITS_EXTENSIONS = [".fits", ".fit", ".fts", ".fz", ".fits.fz"]
IMAGE_EXTENSIONS = [".jpg", ".jpeg", ".png"]


def header_to_wtml_params(header):
    header = header.copy()
    has_cd = "CD1_1" in header
    if has_cd:
        log("\tHeader has CD matrix", level="DEBUG")
    has_crota = "CROTA2" in header
    if has_crota:
        log("\tHeader has CROTA", level="DEBUG")
    has_cdelt = "CDELT1" in header
    if has_cdelt:
        log("\tHeader has CDELT", level="DEBUG")
    has_pc = "PC1_1" in header
    if has_pc:
        log("\tHeader has PC matrix", level="DEBUG")

    wcs = WCS(header)
    crpix = [(header["NAXIS1"] + 1) / 2, (header["NAXIS2"] + 1) / 2]
    crval = wcs.wcs_pix2world(*crpix, 0)
    offset = crpix
    # cd = wh.get_cd(wcs)
    scales, rot, parity = wh.get_scale_rot(header)
    scale = np.sqrt(np.abs(scales[0] * scales[1]))

    return rot, crval, offset, scale, parity


def format_wwt_url(name, rot, crval, offset, scale, bottoms_up, url, thumb_url=''):
    if not isinstance(bottoms_up, bool):
        bottoms_up = bottoms_up > 0
    base_url = "http://www.worldwidetelescope.org/wwtweb/ShowImage.aspx?"
    url_rot = rot + 180
    if url_rot > 180:
        url_rot = -(360 - url_rot)
    options = (
        bottoms_up, 3600 * scale,
        name, url,
        crval[0], crval[1],
        offset[0], offset[1],
        url_rot, thumb_url,
    )
    template = "reverseparity=%s&scale=%.6f&name=%s&imageurl=%s&ra=%.6f&dec=%.6f&x=%.1f&y=%.1f&rotation=%.2f&thumb=%s"
    wwt_url = template  % options
    return base_url + wwt_url

def wwt_url_from_header(header, name, url = None, thumb_url = ""):
    if url is None:
        raise(ValueError("No URL provided"))
    rot, crval, offset, scale, parity = header_to_wtml_params(header)
    bottoms_up = parity > 0
    return format_wwt_url(name, rot, crval, offset, scale, bottoms_up, url)

def create_imageset_dict(name, rot, crval, offset, scale, height, width, url, parity, thumb_url=""):
    log("\n WTML Place & Imageset values: \n", level=3)
    bottoms_up = parity > 0
    full_string = f"""Name: {name} Rotation: {rot:0.2f} Bottoms-up: {bottoms_up} 
                      Center_X: {crval[0]:0.2f} Center_Y: {crval[1]:0.2f} 
                      Offset_X: {offset[0]} Offset_Y: {offset[0]} 
                      Scale: {scale:0.3g} 
                      Height: {height} Width: {width} URL: {url}"""
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
    
    imageset = hc.ImageSet(
        Name = name,
        Url = url,
        Rotation=rot,
        CenterX=crval[0],
        CenterY=crval[1],
        BottomsUp=bottoms_up,
        OffsetX=offset[0],
        OffsetY=offset[1],
    )

    if imageset.Projection == "SkyImage":
        imageset.set_base_degrees_per_tile(scale)
    else:
        imageset.set_base_degrees_per_tile(max(width, height) * scale)
    thumb_url = thumb_url or "https://nova.astrometry.net/image/16942765"
    wwt_url = format_wwt_url(name, rot, crval, offset, scale, bottoms_up, url, )
    
    log("WorldWide Telescope URL:", level=2)
    log(wwt_url, level=2)
    log("\n", level=3)

    return imageset.json()


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


def create_wtml(
    header,
    image_path,
    name="test",
    url=None,
    credits="credits",
    credits_url="credits_url",
    description="description",
    thumbnail_url=None,
    out=None,
    zoom_factor = 1.7,
    place_center = {},
):
    # print out header sayin "Creating WTML file"
    log("\n ****** Creating WTML file ****** \n", level="INFO")
    if out is None:
        out = name + ".wtml"
    log(f"WTML file: {out}", level="INFO")
    
    im = ih.get_PIL_image(image_path)
    width, height = im.size
    if url is None:
        raise ValueError("I need a url")
    
    header = hc.ImageHeader(image_path, header).header
    rot, crval, offset, scale, parity = header_to_wtml_params(header)

    imageset = create_imageset_dict(
        name, rot, crval, offset, scale, height, width, url, parity
    )

    folder = set_folder(name)
    # in WWT zoom is the FOV * 6
    fov = 6 * scale * max(height, width) #in WWT zoom units
    zoom = fov * zoom_factor  # 1.7 is a fudge factor
    # round to three significant digits
    zoom = round(zoom, -int(np.floor(np.log10(abs(zoom)))) + 2)
    place = set_place(name=name, 
                      ra=place_center.get('ra', crval[0]) / 15, 
                      dec=place_center.get('dec', crval[0]), 
                      zoom=place_center.get('zoom', zoom))

    # log('Creating WTML file')
    el_imageset = ET.Element("ImageSet", attrib=imageset)

    el_credits = ET.Element("Credits")
    el_credits.text = credits

    el_credits_url = ET.Element("CreditsUrl")
    el_credits_url.text = credits_url

    el_description = ET.Element("Description")
    el_description.text = description

    el_thumbnail_url = ET.Element("ThumbnailUrl")
    el_thumbnail_url.text = (
        thumbnail_url or image_path
    )

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

    xml.smart_indent_xml(el_folder)
    tree = ET.ElementTree(el_folder)
    tree.write(out)
    xml.split_xml_attributes(out, field="ImageSet")
    return tree, imageset, out

# what do we want to specify when creatng a wtml file?
# name: Place and ImageSet name
# wtml: Name of wtml file
# image_url: where wtml file can find the image
# thumb_url: where wtml file can find the image thumbnail
# outdir: where to save the wtml file
def create_wtml_from_image(
        image_path, 
        wcsfile=None, 
        use_avm = False,
        wtml = None,
        output_dir=None, 
        name=None, 
        image_url=None,
        thumb_url=None,
        suffix="",
        to_github=False,):
    """
    Create a WTML file from an image file.
    image_path: path to image file
    wcsfile: path to wcs file. If not provided it will look for a .wcs file with the same name as the image file
    use_avm: if True, it will try to get the wcs from the AVM tags in the image file
    wtml: name of the wtml file. If not provided it will use the name of the image file
    output_dir: directory where the wtml file will be saved. If not provided it will use the directory of the image file
    name: name of the Place and ImageSet. If not provided it will use the name of the image file
    image_url: url where the wtml file can find the image. If not provided it will use the name of the image file
    thumb_url: url where the wtml file can find the image thumbnail. If not provided it will use the name of the image_url
    suffix: suffix to add to the name of the wtml file (and other outputs)
    """
    log(f"image_path: {image_path}", level="INFO")
    # Get the image
    ext = os.path.splitext(image_path)[1][1:]
    suffix = ih.get_suffix(suffix)

    if name is None:
        name = os.path.splitext(os.path.basename(image_path))[0]
    
    image_header = hc.ImageHeader(image_path, wcsfile=wcsfile, use_avm=use_avm)
    header =image_header.header
    
    # wtml: name of the wtml file
    # wtml_dir: directory where the wtml file will be saved
    if wtml is None:
        wtml = os.path.basename(image_path).replace(f".{ext}", f"{suffix}.wtml")
        wtml_dir = os.path.dirname(image_path)
    else:
        wtml_dir = os.path.dirname(wtml)
        wtml = os.path.basename(wtml)
    
    if output_dir is None:
        if wtml_dir == '':
            output_dir = os.path.dirname(image_path)
        else:
            output_dir = wtml_dir
    
    log(f"output_dir: {os.path.join(output_dir, wtml)}", level="INFO")
    

    if image_url is None:
        image_url = image_path.replace(f".{ext}", f"{suffix}.{ext}")

    log(f"WTML Image URL: {image_url}", level="INFO")
    
    if to_github:
        image_url = ph.to_github(os.path.basename(image_url), web=True)
        
    if thumb_url is None:
        thumb_url = image_url
    
    
    tree, imageset, out = create_wtml(header, image_path, name = name, out=os.path.join(output_dir,wtml), url=image_url, thumbnail_url=thumb_url)
    
    return tree, imageset, out
