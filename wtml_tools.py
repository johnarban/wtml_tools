
from PIL import Image
from astropy.io.fits import Header

from astropy.wcs import WCS
from pyavm import AVM
import os
import numpy as np
from glob import glob
from pyavm.avm import extract_xmp, parse_avm_content
from io import BytesIO
# from lxml.etree import ElementTree as ET
# import lxml.etree as etree
import xml.etree as etree
import xml.etree.ElementTree as ET
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy.wcs.utils import is_proj_plane_distorted, _is_cd_orthogonal
from  aplpy import FITSFigure

import wcs_helpers as wh
from importlib import reload
reload(wh)
from xml_indenter import smart_indent_xml, split_xml_attributes

 

def get_scale_rot(header, force = False):
    # from astrometry.net/net/wcs.py
    wcs = wh.header_to_wcs(header)
    cd = wh.get_cd(wcs = wcs)
    
    scales = wh.get_scale(cd) # in degrees / pixel
    
    rot = wh.get_rot(cd) # in degress
    
    parity = wh.get_parity(cd) # 1 or -1

    print(f'scale: {scales}, rot: {rot:.3f}, parity: {parity}')
    
    return scales, rot, parity


def add_scale_rot(header, from_header = None):
    header = header.copy()
    if from_header is None:
        scale, rotation, parity = get_scale_rot(header)
    else:
        scale, rotation, parity = get_scale_rot(from_header)
    # if parity > 0:
    #     rotation = 180 - rotation
    header['CDELT1'] = scale[0]
    header['CDELT2'] = scale[1]
    header['CROTA2'] = rotation
    # remove CD matrix
    # header['CRPIX2'] = header['NAXIS2'] + 1 - header['CRPIX2']
    return header
 

def do_parity_inversion(header, im):
    # flip image top to bottom 
    # from https://github.com/WorldWideTelescope/wwt-aligner/blob/master/backend/wwt_aligner/driver.py#L449
    parity = wh.get_parity(header = header)
    if parity < 0:
        print('parity inversion not needed')
        # WWT requires images with negative parity
        return header
    print('doing parity inversion')
    do_add_naxisi = 'NAXIS1' in header
    # wcs = WCS(header)
    # hdwork = wcs.to_header(relax=False)
    hdwork = header.copy()
    hdwork['CRPIX2'] = im.height + 1 - hdwork['CRPIX2']
    try:
        hdwork['PC1_2'] *= -1
        hdwork['PC2_2'] *= -1
    except KeyError:
        hdwork['CD1_2'] *= -1
        hdwork['CD2_2'] *= -1
    hdwork = wh.add_NAXES(hdwork,im, add_naxisi = do_add_naxisi)
    return hdwork

# reproject
def get_rot_after_parity_inversion(header):
    header = header.copy()
    header['PC1_2'] *= -1
    header['PC2_2'] *= -1
    return 180 - wh.get_rot(wh.get_cd(header))
    

def rgb_reproject(image, header):
    from_wcs = WCS(header)
    r, g, b  = image.split()
    r = np.array(r)[:,::-1]
    g = np.array(g)[:,::-1]
    b = np.array(b)[:,::-1]
    to_wcs,shape_out = find_optimal_celestial_wcs([(r.shape,from_wcs)])
    r, footprint = reproject_interp((r, from_wcs), to_wcs, shape_out = shape_out)
    g, footprint = reproject_interp((g, from_wcs), to_wcs, shape_out = shape_out)
    b, footprint = reproject_interp((b, from_wcs), to_wcs, shape_out = shape_out)
    out = np.nan_to_num([r,g,b])
    return out, to_wcs
    
    
   
   
def output_avm(im, header, name='test', suffix = '', path_out='./', ext='jpg',  add_rot = False):
    
    if len(suffix) > 0 and '_' != suffix[0]:
        suffix = '_' + suffix
    name = name + suffix
    
    path_out = path_out.rstrip('/')
    out_tagged = path_out + '/' + name + '_tagged.' + ext
    print(f'outputting {name} to {out_tagged}')
    
    # make sure we a have a valid JPEG for AVM to work with
    temp = 'temp_can_delete.' + ext # temporary_file
    im.save(temp)
    
    
    # hdr = wcs.to_header()
    # hdr = wh.add_NAXES(hdr,im)
    
    # if add_rot:
    #     hdr = add_scale_rot(hdr)
    #     wcs = WCS(hdr, relax = True)
    
    if isinstance(header, WCS):
        avm = AVM.from_wcs(header, shape = (im.height, im.width))
    else:
        avm = AVM.from_header(header)
    
    avm.embed(temp,out_tagged)
    
    return out_tagged, avm
    


def clean_image(filename):
    im = Image.open(filename)
    im.save(filename)
    # # remove first line of a file
    # with open(filename, 'r') as fin:
    #     data = fin.read().splitlines(True)
    # with open(filename, 'w') as fout:
    #     fout.writelines(data[1:])

def get_header(wcsfile, clean = False, clean_sip = False, clean_all = False):
    header = Header.fromfile(wcsfile)
    if clean or clean_all:
        header = wh.clean_header(header)
    if clean_sip or clean_all:
        header = wh.remove_sip(header)
    return header

def add_avm_tags_simple(image, name = None, wcsfile = None, use_inversion = True, add_rot = False, suffix = '', path_out = '.'):
    im = Image.open(image)
    ext = image.split('.')[-1]
    
    if len(suffix) > 0 and '_' != suffix[0]:
        suffix = '_' + suffix
    
    if wcsfile is None:
        wcsfile = image.replace('.'+ext,'.wcs')
    
    if name is None:
        name = image.replace(f'.{ext}','')

    header = get_header(wcsfile, clean_all=True)
    header = wh.add_NAXES(header,im, add_naxisi = True)
    print('Original parity: ', wh.get_parity(header = header))

    wcs = WCS(header)
    
    if use_inversion:
        new_header = do_parity_inversion(header,im) 
        if add_rot:
            rot = get_rot_after_parity_inversion(header)
            parity = wh.get_parity(header)
            new_header['CROTA2'] = parity * rot
            scale = wh.get_scale(new_header)
            new_header['CDELT1'] = scale[0]
            new_header['CDELT2'] = scale[1]
        wcs_for_embed = WCS(new_header)
        print('New parity: ', wh.get_parity(header = header))

    else:
        new_header = add_scale_rot(header)
        new_header = wh.remove_cd(new_header)
        wcs_for_embed = WCS(new_header)
    out_tagged, avm = output_avm(im, wcs_for_embed, 
                            name=name,
                            suffix = suffix, 
                            path_out=path_out, 
                            add_rot = add_rot,
                            ext=ext)
    return header, wcs, new_header, wcs_for_embed, avm, out_tagged, im

def add_avm_reproject(image, wcsfile = None, path_out = ''):
    im = Image.open(image)
    if wcsfile is None:
        wcsfile = image.replace('jpg','wcs')

    header = wh.clean_header(Header.fromfile(wcsfile))
    header = wh.add_NAXES(header,im, add_naxisi = True)
    
    # reproject
    new_im, new_wcs = rgb_reproject(im, header)
    new_image = Image.fromarray(new_im.astype('uint8').transpose(1,2,0))
    out_tagged, avm = output_avm(new_image, new_wcs, 
                            name=image.replace('.jpg','_reproj'), 
                            ext='jpg', path_out=path_out)
    return header, WCS(header), new_wcs, avm, out_tagged, im, new_image
    
    
    
    
def embed_avm_from_wcs(image, name = None, wcsfile = None, path_out = '', use_inversion = True, suffix = '', add_rot = False):
    im = Image.open(image)
    if wcsfile is None:
        wcsfile = image.replace('jpg','wcs')
    if name is None:
        name = image.replace('.jpg','')
    
    header = Header.fromfile(wcsfile)
    header = wh.clean_header(header)
    header = wh.remove_sip(header)
    header = wh.add_NAXES(header,im, add_naxisi = True)
    # header['CRPIX2'] = im.height + 1 - header['CRPIX2']
    wcs = WCS(header)
    
    out_tagged, avm = output_avm(im, wcs, 
                            name=name,
                            suffix = suffix, 
                            ext='jpg', path_out=path_out, add_rot = add_rot)
    return header, wcs, avm, out_tagged, im


def preview_image(image,coords=None):
    g =FITSFigure(image)
    g.add_grid()
    g.set_theme('pretty')
    if coords is not None:
        g.show_markers(coords.ra,coords.dec,marker='x',facecolor='red',edgecolor='red',s=15)
    g.axis_labels.show()
    # g.tick_labels.set_xformat('hh:mm:ss.ss')
    g.tick_labels.show()
    g.tick_labels.set_style('colons')
    g.show_rgb()
    return g
    
    

# print out root and children as a tree showing attributes and values, recursively
def print_tree(node, level=0, print_attributes=False):
    # print('-'*20)
    print('  '*level + node.tag)
    # print the attributes
    if print_attributes:
        for key, value in node.attrib.items():
            print('  '*(level+1) + key,':', value)
    for child in node:
        print_tree(child, level+1, print_attributes=print_attributes)
        

def create_imageset_dict(name, rot, crval, offset, scale, height, width, url, parity, cd):
    print('  create_imageset_dict')
    print('  name',name)
    print('  parity',parity)
    bottoms_up = parity < 0
    if bottoms_up:
        # in negative parity the rotation needs to be reverse 
        rot = -rot   
    print('  rot',rot)
    print('  crval',crval)
    print('  offset',offset)
    print('  scale',scale)
    print('  height',height)
    print('  width',width)
    # repository is the github repo in the url
    repository = '/'.join(url.split('/')[3:5])
    print(f'  repository @{repository}')
    directory = '/'.join(url.split('/')[4:-1])
    print(f'  directory: {directory}/')
    file = url.split('/')[-1]
    print(f'  image: {file}')
    
    print(f'  url: {url}\n')
     

    
    imageset = {
        "DemUrl":"",
        "MSRCommunityId":"0",
        "MSRComponentId":"0",
        "Permission":"0",
        'Generic': 'False',
        'DataSetType': 'Sky',
        'BandPass': 'Visible',
        'Url': url,
        'TileLevels': '0',
        'WidthFactor': '2',
        'Sparse': 'True',
        'Rotation': str(rot),
        'QuadTreeMap': '',
        'Projection': 'SkyImage',
        'Name': name,
        'FileType': '.jpg',
        'CenterY': str(crval[1]), # ra in decimal degrees
        'CenterX': str(crval[0]), # dec in decimal degrees
        'BottomsUp':str(bottoms_up), # True if it should be inverted
        'StockSet': 'False',
        'ElevationModel': 'False',
        'OffsetX': str((offset[0])),  # ra offset in pixels
        'OffsetY': str((offset[1])), # dec offset in pixels
        'BaseTileLevel': '0',
        'BaseDegreesPerTile': '',
        'ReferenceFrame': '',
        'MeanRadius': '0'
    }
    
    if imageset['Projection'] == 'SkyImage':
        imageset['BaseDegreesPerTile'] =  str(np.sqrt(cd[0,1]**2 + cd[1,1]**2)) # scal_y
    else:
        imageset['BaseDegreesPerTile'] =  str(max(width,height) * scale) #str(np.sqrt(cd[0,1]**2 + cd[1,1]**2)) #str(height_arcmin * max(height,width) / (60 * height))
    
    base_url = "http://www.worldwidetelescope.org/wwtweb/ShowImage.aspx?"
    options = (bottoms_up, 3600 * scale, name, url, crval[0], crval[1], offset[0], offset[1], rot+180, 'https://nova.astrometry.net/image/16942765')
    wwt_url = 'reverseparity=%s&scale=%.6f&name=%s&imageurl=%s&ra=%.6f&dec=%.6f&x=%.1f&y=%.1f&rotation=%.2f&thumb=%s' % options
    print('WorldWide Telescope URL:')
    print(base_url + wwt_url)
    print('\n')

    return imageset

def set_folder(name):
    folder = {
        "MSRCommunityId":"0",
        "MSRComponentId":"0",
        "Permission":"0",
        "Name":name,
        "Group":"Explorer",
        "Searchable":"False",
    }
    return folder 
    
def set_place(name, ra, dec, rotation = 0, zoom = 10):
    place = {
        "Name":name,
        "DataSetType":"Sky",
        "RA": str(ra), # ra in hours
        "Dec": str(dec), # dec in decimal degrees
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
        
def set_xml_element(el, attribute = None, text = None):
    if attribute is not None:
        try:
            el.attrib = attribute
        except:
            for key, value in attribute.items():
                el.set(key, value)
            
    if text is not None:
        el.text = text
    
def create_wtml(header,im,
                name='test',
                url = None,
                credits = 'credits', 
                credits_url = 'credits_url',
                thumbnail_url = None,
                out=None):
    # print out header sayin "Creating WTML file"
    print('\n ****** Creating WTML file ****** \n')
    if out is None:
        out = name + '.wtml'
    print('writing to',out)
    wcs = WCS(header)
    # move reference to center of the image
    # get the specs we need for the wtml file
    # crpix = header['CRPIX1'],header['CRPIX2']
    # crval = header['CRVAL1'],header['CRVAL2']
    crpix = [(header['NAXIS1']+1)/2,(header['NAXIS2']+1)/2]
    crval =  wcs.wcs_pix2world(*crpix,1)
    offset = crpix
    # crpix = [0,0]
    # crval = wcs.wcs_pix2world([[crpix[0],crpix[1]]],1,)[0]
    # offset =  0,(header['NAXIS2']+1)
    
    cd = wh.get_cd(wcs)
    scales, rot, parity = get_scale_rot(header)
    scale = np.sqrt(np.abs(scales[0]*scales[1]))
    if url is None:
        print('I need a url')
        return 0
    url = url
    imageset = create_imageset_dict(name,
                                        rot,
                                        crval,
                                        offset,
                                        scale,
                                        im.height,
                                        im.width,
                                        url,
                                        parity,
                                        cd)
    folder = set_folder(name)
    ra, dec = wcs.wcs_pix2world((header['NAXIS1']-1)/2,(header['NAXIS2']-1)/2,1)
    place = set_place(name,ra/15, dec)
    # add imageset to xml
    # tree = parse_xml('blank.wtml')
    # root = tree.getroot()
    el_imageset = ET.Element('ImageSet', attrib = imageset)
    el_credits = ET.Element('Credits')
    el_credits.text = credits
    el_credits_url = ET.Element('CreditsUrl')
    el_credits_url.text = credits_url
    el_description = ET.Element('Description')
    el_description.text = 'description'
    el_thumbnail_url = ET.Element('ThumbnailUrl')
    el_thumbnail_url.text = thumbnail_url or 'https://nova.astrometry.net/image/16942765'
    
    
    el_imageset.append(el_description)
    el_imageset.append(el_credits)
    el_imageset.append(el_credits_url)
    el_imageset.append(el_thumbnail_url)
    
    el_forground = ET.Element('ForegroundImageSet')
    el_forground.append(el_imageset)
    
    el_place = ET.Element('Place', attrib = place)
    el_place.append(el_forground)
    
    el_folder = ET.Element('Folder', attrib = folder)
    el_folder.append(el_place)
    
    smart_indent_xml(el_folder)
    tree = ET.ElementTree(el_folder)
    tree.write(out)
    split_xml_attributes(out,field = 'ImageSet')
    return tree, imageset