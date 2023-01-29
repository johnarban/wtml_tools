from astropy.io.fits import Header
from astropy.wcs import WCS
from math import atan2, atan
from astropy.wcs.utils import _is_cd_orthogonal
from numpy import pi


 # The functions defined in the file are as follows:
# - get_cd
# - get_parity
# - header_to_wcs
# - get_quadrant
# - get_rot
# - get_scale
# - clean_header
# - remove_sip
# - add_NAXES
# - remove_cd
# - header_cd_to_cdelt_crota

# HEADER/WCS manipulation
def get_cd(header = None, wcs = None):
    """Return the CD matrix from a header or wcs object.
    This is a convenience function that will return the CD matrix
    from either a header or a wcs object.
    This works using the astropy.wcs.WCS.pixel_scale_matrix property.
    This property derives from either wcs.wcs.cd or wcs.wcs.pc & wcs.wcs.cdelt.
    """
    # check if user was lazy and header is actually a wcs object
    if (isinstance(header, WCS)):
        wcs = header
        header = None
    
    if wcs is not None:
        return wcs.pixel_scale_matrix
    else:
        return WCS(header).pixel_scale_matrix


def det_cd(cd_matrix):
    """Return the determinant of a CD matrix.
    
    the cd Matrix looks like [ CD1_1 CD1_2 ]
                             [ CD2_1 CD2_2 ]
    written in the form of a numpy array:
                             [ [0,0] [0,1] ]
                             [ [1,0] [1,1] ]
    """
    return cd_matrix[0,0] * cd_matrix[1,1] - cd_matrix[0,1] * cd_matrix[1,0]



def get_parity(cd = None, header = None):
    """get's the parity of the cd matrix
    FITS typically have positive parity
    while JPEGS tend to have negative parity
    See https://github.com/WorldWideTelescope/toasty/blob/master/toasty/builder.py#L290
    for a discussion on parity
    """
    if cd is None:
        cd = get_cd(header)

    if (det_cd(cd) >= 0) or _is_cd_orthogonal(cd,1e-5):
        return  1
    else:
        return -1

def header_to_wcs(supposed_wcs):
    """Converts a header to a wcs object if it isn't already one.
    """
    if isinstance(supposed_wcs, WCS):
        return supposed_wcs
    else:
        return WCS(supposed_wcs)


def get_quadrant(x, y):
    pos_x = x > 0
    pos_y = y > 0
    if pos_x and pos_y:
        return 'I'
    elif not pos_x and pos_y:
        return 'II'
    elif not pos_x and not pos_y:
        return 'III'
    elif pos_x and not pos_y:
        return 'IV'
    else:
        raise ValueError('x and y are zero!')


def get_rot(cd, wwt = False):
    """Return the rotation angle of the CD matrix.
    This can be measured two ways:
    === Method 1 === (@astrometry.net/net/wcs.py#L80)
      T = parity * cd[0,0] + cd[1,1]
      A = parity * cd[1,0] - cd[0,1]
      rot = -arctan2( A, T) * 180 / pi
    
    === Method 2 === (@wwt_data_formats/wwt_data_formats/imageset.py#L479)
        T =  cd[1,1]
        A = cd[0,1]
        rot = arctan2( - parity A, - T) * 180 / pi
    
    We use method 2 because it is what WWT uses.
    """
    parity = get_parity(cd)
    log = lambda A, T: print(f'angle {atan(A/T)*180/pi:0.3f} in quadrant {get_quadrant(T,A)}')
    if wwt:
        T = cd[1,1]
        A = cd[0,1]
        log(A,T)
        return atan2(- parity * A, -T) * 180 / pi
    else:
        parity = get_parity(cd)
        T = parity * cd[0,0] + cd[1,1]
        A = parity * cd[1,0] - cd[0,1]
        log(A,T)
        return -atan2(- parity * A, -T) * 180 / pi
        
def get_scale(cd):
    """Return the scale of the CD matrix.
    """
    return (cd**2).sum(axis=0)**0.5 

    

def clean_header(header, inplace = False, verbose = False):
    """Remove all HISTORY and COMMENT cards from a header.
    """
    if verbose: print('cleaning header')
    if not inplace:
        header = header.copy()
    header.remove('HISTORY',remove_all=True,ignore_missing=True)
    header.remove('COMMENT',remove_all=True,ignore_missing=True)
    return header

def remove_sip(header,inplace = False, verbose = False):
    """Remove all SIP related keywords from a header.
    """
    if verbose: print('removing SIP')
    if not inplace:
        header = header.copy()
    header['CTYPE1'] = header['CTYPE1'].replace('-SIP','')
    header['CTYPE2'] = header['CTYPE2'].replace('-SIP','')
    for i in [0,1,2]:
        for j in [0,1,2]:
            header.remove(f'A_{i}_{j}',ignore_missing=True)
            header.remove(f'B_{i}_{j}',ignore_missing=True)
            header.remove(f'AP_{i}_{j}',ignore_missing=True)
            header.remove(f'BP_{i}_{j}',ignore_missing=True)
    i = 'ORDER'
    header.remove(f'A_{i}',ignore_missing=True)
    header.remove(f'B_{i}',ignore_missing=True)
    header.remove(f'AP_{i}',ignore_missing=True)
    header.remove(f'BP_{i}',ignore_missing=True)
    return header

# fixup headers
def add_NAXES(header,im, add_naxisi = True, inplace = False, verbose = False):
    """Add NAXIS and NAXIS1, NAXIS2 keywords to a header.
    """
    if verbose: print('adding NAXES')
    if not inplace:
        header = header.copy()
    header['NAXIS'] = 2
    if add_naxisi:
        header['NAXIS1'] = im.width
        header['NAXIS2'] = im.height
    return header


def remove_cd(header, verbose = False):
    if verbose: print('removing CD matrix')
    header = header.copy()
    if hasattr(WCS(header).wcs,'cd'):
        pre = 'CD'
    elif hasattr(WCS(header).wcs,'pc'):
        pre = 'PC'
    else:
        return header
    header.remove(pre + '1_1',ignore_missing=True)
    header.remove(pre + '1_2',ignore_missing=True)
    header.remove(pre + '2_1',ignore_missing=True)
    header.remove(pre + '2_2',ignore_missing=True)
    return header
    
def header_cd_to_cdelt_crota(header):
    """
    Convert a header with a CD matrix to a header with CDELT and CROTA keywords.
    """
    header = header.copy()
    cd = get_cd(header = header)
    header = remove_cd(header)
    
    parity = get_parity(header = header)
    
    cdelt = get_scale(cd)
    crota = get_rot(cd)
    if parity > 0:
        # just need parity inversion
        cdelt[1] *= -1
        crota = crota
    header['CDELT1'] = cdelt[0]
    header['CDELT2'] = cdelt[1]
    header['CROTA2'] = crota 
    return header
    
