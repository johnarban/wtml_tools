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

import logger as logger
reload(logger)
FITS_EXTENSIONS = [".fits", ".fit", ".fts", ".fz", ".fits.fz"]
IMAGE_EXTENSIONS = [".jpg", ".jpeg", ".png"]


def header_to_wtml_params(header, for_url = False, nudge={'x':0, 'y':0}):
    header = header.copy()
    has_cd = "CD1_1" in header
    if has_cd:
        logger.log("\tHeader has CD matrix", level="DEBUG")
    has_crota = "CROTA2" in header
    if has_crota:
        logger.log("\tHeader has CROTA", level="DEBUG")
    has_cdelt = "CDELT1" in header
    if has_cdelt:
        logger.log("\tHeader has CDELT", level="DEBUG")
    has_pc = "PC1_1" in header
    if has_pc:
        logger.log("\tHeader has PC matrix", level="DEBUG")

    wcs = WCS(header)
    crpix = [(header["NAXIS1"] + 1) / 2, (header["NAXIS2"] + 1) / 2]
    crval = wcs.wcs_pix2world(*crpix, 1)
    if for_url:
        offset = crpix[0] + nudge['x'], crpix[1] + nudge['y']
    else:
        offset = crpix[0] + nudge['x'], crpix[1] + nudge['y']
    
    
    # cd = wh.get_cd(wcs)
    scales, rot, parity = wh.get_scale_rot(header)
    scale = np.sqrt(np.abs(scales[0] * scales[1]))

    return {'rot':rot, 'crval':crval, 'offset':offset, 'scale':scale, 'parity':parity}


def format_wwt_url(name, rot = None, crval = None, offset = None, scale = None, reverse_parity = None, url = '', thumb_url='', **kwargs):
    if not isinstance(reverse_parity, bool):
        reverse_parity = reverse_parity > 0
    base_url = "http://www.worldwidetelescope.org/wwtweb/ShowImage.aspx?"
    url_rot = rot + 180 # https://wwt-api-client.readthedocs.io/en/latest/endpoints/legacy/show-image.html#endpoint-showimage
    if url_rot > 360:
        url_rot = url_rot - 360
    options = (
        reverse_parity, 3600 * scale,
        name, url,
        crval[0], crval[1],
        offset[0], offset[1],
        url_rot, thumb_url,
    )
    logger.log('\n ****** Creating WTML URL ****** \n', level="INFO")
    logger.log(f"\tReverse Parity: {reverse_parity}", level="STATS")
    logger.log(f"\tScale: {3600 * scale} arcsec", level="STATS")
    logger.log(f"\tName: {name}", level="STATS")
    logger.log(f"\tURL: {url}", level="STATS")
    logger.log(f"\tRA: {crval[0]}", level="STATS")
    logger.log(f"\tDec: {crval[1]}", level="STATS")
    logger.log(f"\tX: {offset[0]}", level="STATS")
    logger.log(f"\tY: {offset[1]}", level="STATS")
    logger.log(f"\tRotation: {url_rot}", level="STATS")
    template = "reverseparity=%s&scale=%.6f&name=%s&imageurl=%s&ra=%.6f&dec=%.6f&x=%.1f&y=%.1f&rotation=%.2f&thumb=%s"
    wwt_url = template  % options
    logger.log(f"\nWorldWide Telescope URL:\n\t {base_url + wwt_url}\n", level="INFO")
    return base_url + wwt_url

def wwt_url_from_header(header, name, url = None, thumb_url = "", nudge = {'x':0, 'y':0}):
    if url is None:
        raise(ValueError("No URL provided"))
    wtml_params = header_to_wtml_params(header, for_url=True, nudge = nudge)
    reverse_parity = wtml_params['parity'] > 0
    return format_wwt_url(name, **wtml_params, reverse_parity = reverse_parity, url = url)

def create_imageset_dict(name = '', rot = None, crval = None, offset = None, scale = None, height = None, width = None, url = None, parity = None, thumb_url=""):
    logger.log("\n WTML Place & Imageset values: \n", level="STATS")
    bottoms_up = parity > 0
    full_string = f"""Name: {name} 
                      Rotation: {rot:0.2f}
                      Bottoms-up: {bottoms_up} 
                      Center_X: {crval[0]:0.7f} Center_Y: {crval[1]:0.7f} 
                      Offset_X: {offset[0]} Offset_Y: {offset[1]} 
                      Scale: {scale:0.3g} Height: {height} Width: {width} 
                      URL: {url}"""
    logger.log(full_string, level="STATS")
    if "github" in url:
        # repository is the github repo in the url
        repository = "/".join(url.split("/")[3:5])
        logger.log(f"  repository @{repository}", level="DEBUG")
        directory = "/".join(url.split("/")[4:-1])
        logger.log(f"  directory: {directory}/", level="DEBUG")
    file = url.split("/")[-1]
    logger.log(f"\timage: {file}", level="DEBUG")

    logger.log(f"\turl: {url}\n", level="DEBUG")
    
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
    nudge = {'x':0, 'y':0},
):
    # print out header sayin "Creating WTML file"
    logger.log("\n ****** Creating WTML file ****** \n", level="INFO")
    if out is None:
        out = name + ".wtml"
    logger.log(f"WTML file: {out}", level="INFO")
    
    im = ih.get_PIL_image(image_path)
    width, height = im.size
    if url is None:
        raise ValueError("I need a url")
    
    header = hc.ImageHeader(image_path, header).header
    wtml_params = header_to_wtml_params(header, nudge=nudge)

    imageset = create_imageset_dict(name, **wtml_params, height = height, width = width, url = url)
    wwt_url = wwt_url_from_header(header, name, url, thumb_url=thumbnail_url or image_path, nudge=nudge)
    
    folder = set_folder(name)
    # in WWT zoom is the FOV * 6
    fov = 6 * wtml_params['scale'] * max(height, width) #in WWT zoom units
    zoom = fov * zoom_factor  # 1.7 is a fudge factor
    # round to three significant digits
    zoom = round(zoom, -int(np.floor(np.log10(abs(zoom)))) + 2)
    place = set_place(name=name, 
                      ra=place_center.get('ra', wtml_params['crval'][0]) / 15, 
                      dec=place_center.get('dec', wtml_params['crval'][1]), 
                      zoom=place_center.get('zoom', zoom))

    # logger.log('Creating WTML file')
    el_imageset = ET.Element("ImageSet", attrib=imageset)

    el_credits = ET.Element("Credits")
    el_credits.text = credits

    el_credits_url = ET.Element("CreditsUrl")
    el_credits_url.text = credits_url
    
    try:
        el_description = ET.Element("Description")
        # interpret description as raw html
        el_description.text = ET.CDATA(description)
    except:
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

def add_description_credits_to_wtml(wtml_file, 
                                    description = "description", 
                                    credits = "credits", 
                                    credits_url = "credits_url",
                                    place_name = None,
                                    imageset_name = None
                                    ):
    # wtml is an xml file
    # it shoulud have <Description> and <Credits> and <CreditsUrl> tags
    # we want to add the description and credits to the wtml file
    
    # parse the xml
    try:
        tree = ET.parse(wtml_file)
    except:
        print('error',wtml_file)
        return
    root = tree.getroot()
    
    imageset = root.find(".//ImageSet")
    
    append_description = False
    append_credits = False
    append_url = False
    
    # find the <Description> tag
    description_tag = root.find(".//Description")
    if description_tag is None:
        append_description = True
        # if there is no description tag, create one
        description_tag = ET.Element("Description")
    # add the description to the tag
    description_tag.text = description
    
    # find the <Credits> tag
    credits_tag = root.find(".//Credits")
    if credits_tag is None:
        append_credits = True
        # if there is no credits tag, create one
        credits_tag = ET.Element("Credits")
        imageset.append(credits_tag)
    # add the credits to the tag
    credits_tag.text = credits
    
    # find the <CreditsUrl> tag
    credits_url_tag = root.find(".//CreditsUrl")
    if credits_url_tag is None:
        append_url = True
        # if there is no credits_url tag, create one
        credits_url_tag = ET.Element("CreditsUrl")
        imageset.append(credits_url_tag)
    # add the credits_url to the tag
    credits_url_tag.text = credits_url
    
    
    if place_name  is not None:
        place = root.find(".//Place")
        place.attrib['Name'] = place_name
    
    if imageset_name is not None:
        imageset.attrib['Name'] = imageset_name
        
    if append_description:
        imageset.append(description_tag)
    if append_credits:
        imageset.append(credits_tag)
    if append_url:
        imageset.append(credits_url_tag)
    
    #print the xml file
    xml.smart_indent_xml(root)
    # write the xml to the wtml file
    tree.write(wtml_file)
    xml.split_xml_attributes(wtml_file, field="ImageSet")
    # replace the &lt; and &gt; with < and >
    xml.replace_lt_gt(wtml_file)
    
    # repla

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
        to_github=False,
        write_avm=False,
        write_header=False,
        zoom_factor = 1.7,
        nudge={'x':0, 'y':0},
        rotate_cd_matrix_by = 0,
        force_avm_180 = False,
        force_image_center = False,
        force_bottoms_up = False,
        to_fits=False,
        mods = None,
        place_center= {},
        description = "description",
        credits = "credits",
        credits_url = "credits_url",
        ):
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
    foce_avm_180: if True, it will force the AVM rotatiion to be angle - 180 degrees. 
    """
    logger.log(f"image_path: {image_path}", level="INFO")
    # Get the image
    ext = os.path.splitext(image_path)[1][1:]
    suffix = ih.get_suffix(suffix)

    if name is None:
        name = os.path.splitext(os.path.basename(image_path))[0]
    
    if wcsfile is None:
        logger.log(f"wcsfile not provided.", level="INFO")
    else:
        logger.log(f"wcsfile: {wcsfile}", level="INFO")
    image_header = hc.ImageHeader(image_path, 
                                  wcsfile=wcsfile, 
                                  use_avm=use_avm, 
                                  rotate_cd_matrix = rotate_cd_matrix_by,
                                  force_image_center = force_image_center,
    )
    if force_bottoms_up:
        image_header.header = wh.flip_parity(image_header.header)
    
    if mods is not None:
        image_header.header = wh.modify_header(image_header.header, **mods)
    
    header = image_header.header
    
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
    
    logger.log(f"output_dir: {os.path.join(output_dir, wtml)}", level="INFO")
    

    if image_url is None:
        image_url = image_path.replace(f".{ext}", f"{suffix}.{ext}")

    logger.log(f"WTML Image URL: {image_url}", level="INFO")
    
    if to_github:
        image_url = ph.to_github(os.path.basename(image_url), web=True)
        
    if thumb_url is None:
        thumb_url = image_url
    
    if write_avm:
        logger.log(f"Writing AVM to {image_path}", level="INFO")
        image_header.write_avm(output_dir, force_180=False)
    
    # save the header
    if write_header:
        logger.log(f"Writing header to {image_path}", level="INFO")
        image_header.write_header(name = 'header.hdr')
    
    if to_fits:
        logger.log(f"Writing FITS to {image_path}", level="INFO")
        image_header.write_fits(output_dir, name='image.fits')
        
    tree, imageset, out = create_wtml(header, 
                                      image_path, 
                                      name = name, 
                                      out=os.path.join(output_dir,wtml), 
                                      url=image_url, 
                                      thumbnail_url=thumb_url, 
                                      zoom_factor=zoom_factor, 
                                      nudge=nudge, 
                                      place_center=place_center,
                                      description=description,
                                      credits=credits,
                                      credits_url=credits_url,
                                      )
    
    return tree, imageset, out, image_header
