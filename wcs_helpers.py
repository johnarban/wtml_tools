from astropy.io.fits import Header
from astropy.wcs import WCS
from math import atan2, atan
from astropy.wcs.utils import _is_cd_orthogonal
from numpy import pi

DEBUG_LEVELS = {"DEBUG": 0, "INFO": 1}
DEBUG_LEVEL = 1


def set_debug_level(level):
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


# The functions defined in the file are as follows:
# - get_cd
# - get_cd_sign
# - header_to_wcs
# - get_quadrant
# - get_rot
# - get_scale
# - clean_header
# - remove_sip
# - add_NAXES
# - remove_cd
# - header_cd_to_cdelt_crota


## CD Matrix tools
## FITS originally used CDELT and CROTA to describe the scale and rotation
## Then it switched to using
def get_cd(header=None, wcs=None):
    """Return the CD matrix from a header or wcs object.
    This is a convenience function that will return the CD matrix
    from either a header or a wcs object.
    This works using the astropy.wcs.WCS.pixel_scale_matrix property.
    This property derives from either wcs.wcs.cd or wcs.wcs.pc & wcs.wcs.cdelt.
    """
    # check if user was lazy and header is actually a wcs object
    if isinstance(header, WCS):
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
    return cd_matrix[0, 0] * cd_matrix[1, 1] - cd_matrix[0, 1] * cd_matrix[1, 0]


def pretty_print_matrix(matrix, name=""):
    """Prints a matrix in a pretty way."""
    str1 = f" [ {matrix[0,0]:9.3g} {matrix[0,1]:9.3g} ]"
    str2 = f" [ {matrix[1,0]:9.3g} {matrix[1,1]:9.3g} ]"
    return str1 + "\n" + str2


def get_parity(cd=None, header=None):
    """Return the parity of the CD matrix."""
    return -get_cd_sign(cd=cd, header=header)


def is_JPEGLike(parity=None, cd_sign=None):
    """Return True if the parity is negative."""
    if cd_sign is not None:
        return cd_sign > 0
    return parity < 0


def get_cd_sign(cd=None, header=None):
    """get's the cd_sign of the cd matrix
    FITS typically have positive cd_sign
    while JPEGS tend to have negative cd_sign
    See https://github.com/WorldWideTelescope/toasty/blob/master/toasty/builder.py#L290
    for a discussion on cd_sign
    """
    if cd is None:
        cd = get_cd(header)

    if det_cd(cd) >= 0:
        return 1
    else:
        return -1


def flip_parity(header, height=None, inplace=False):
    """
    Flip the parity of the FITS header
    This cannot be done with a header that has only
    scale and rotation matrix, it must have a CD or PC matrix.
    We flip the parity by flipping the 1_2 and 2_2 elements of
    PC/CD matrix
    the height should be from NAXIS2, or given
    """
    if not inplace:
        header = header.copy()
    original_cd_sign = get_cd_sign(header=header)
    original_pixel_scale_matrix = get_cd(header)
    if height is None:
        if "NAXIS2" in header:
            height = header["NAXIS2"]
        else:
            raise ValueError("Cannot flip parity without height")

    header["CRPIX2"] = height + 1 - header["CRPIX2"]
    if "PC1_2" in header:
        log("flipping PC")
        header["PC1_2"] *= -1
        header["PC2_2"] *= -1
    elif "CD1_2" in header:
        log("flipping CD")
        header["CD1_2"] *= -1
        header["CD2_2"] *= -1
    else:
        raise ValueError("Cannot flip parity without CD or PC matrix")
    new_cd_sign = get_cd_sign(header=header)
    new_pixel_scale_matrix = get_cd(header)

    # pretty print the matrices showing the change
    log(f"Original CD matrix: cd_sign = {original_cd_sign}", level=1)
    log(pretty_print_matrix(original_pixel_scale_matrix), level=0)
    log(f"New CD matrix: cd_sign = {new_cd_sign}", level=1)
    log(pretty_print_matrix(new_pixel_scale_matrix), level=0)

    return header


def header_to_wcs(supposed_wcs):
    """Converts a header to a wcs object if it isn't already one."""
    if isinstance(supposed_wcs, WCS):
        return supposed_wcs
    else:
        return WCS(supposed_wcs)


def get_quadrant(x, y):
    pos_x = x > 0
    pos_y = y > 0
    if pos_x and pos_y:
        return "I"
    elif not pos_x and pos_y:
        return "II"
    elif not pos_x and not pos_y:
        return "III"
    elif pos_x and not pos_y:
        return "IV"
    else:
        raise ValueError("x and y are zero!")


def positive_angle(angle, radians=False):
    if radians:
        return (angle + 2 * pi) % (2 * pi)
    else:
        return (angle + 360) % 360


def get_rot_from_cd(cd, wwt=False):
    """Return the rotation angle of the CD matrix.
    This can be measured two ways:
    === Method 1 (preferred) === (@astrometry.net/net/wcs.py#L80)
      T = cd_sign * cd[0,0] + cd[1,1]
      A = cd_sign * cd[1,0] - cd[0,1]
      rot = -arctan2( A, T) * 180 / pi

    === Method 2 === (@wwt_data_formats/wwt_data_formats/imageset.py#L479)
        T =  cd[1,1]
        A = cd[0,1]
        rot = arctan2( - cd_sign A, - T) * 180 / pi

    We use method 2 because it is what WWT uses.
    """
    log("Finding the rotation angle")
    cd_sign = get_cd_sign(cd)
    log(f"\tCD sign ={cd_sign}", level = 'DEBUG')
    # wwt
    T = cd[1, 1]
    A = cd[0, 1]
    rot_wwt = atan2(-cd_sign * A, -T) * 180 / pi
    log(f"\tRot (WWT) = {rot_wwt:0.3f} degrees", level = 'DEBUG')

    # astrometry.net
    T = cd_sign * cd[0, 0] + cd[1, 1]
    A = cd_sign * cd[1, 0] - cd[0, 1]
    rot_astrom = -atan2(-cd_sign * A, -T) * 180 / pi
    log(f"\tRot (Astrometry.net) = {rot_astrom:0.3f} degrees", level = 'DEBUG')
    if wwt:
        log(f"\tUsing WWT method. Angle :{rot_wwt:0.3f}", level='INFO')
        return rot_wwt
    else:
        log(f"\tUsing Astrometry.net method. Angle :{rot_astrom:0.3f}", level='INFO')
        return rot_astrom


def get_scale(cd):
    """Return the scale of the CD matrix."""
    return (cd**2).sum(axis=0) ** 0.5


def get_scale_rot(header, force=False, wwt=False):
    # from astrometry.net/net/wcs.py
    wcs = header_to_wcs(header)
    cd = get_cd(wcs=wcs)

    scales = get_scale(cd)  # in degrees / pixel

    rot = get_rot_from_cd(cd, wwt=wwt)  # in degress

    parity = get_parity(cd)  # 1 or -1

    log("get_scale_rot::")
    log(f"\tscale: {scales}\n\trot: {rot:.3f}\n\tparity: {parity}", level=0)

    return scales, rot, parity


def clean_header(header, inplace=False, verbose=False):
    """Remove all HISTORY and COMMENT cards from a header."""
    log("cleaning header", level=0)
    if not inplace:
        header = header.copy()
    header.remove("HISTORY", remove_all=True, ignore_missing=True)
    header.remove("COMMENT", remove_all=True, ignore_missing=True)
    return header


def remove_sip(header, inplace=False, verbose=False):
    """Remove all SIP related keywords from a header."""
    if verbose:
        log("removing SIP", level=0)
    if not inplace:
        header = header.copy()
    header["CTYPE1"] = header["CTYPE1"].replace("-SIP", "")
    header["CTYPE2"] = header["CTYPE2"].replace("-SIP", "")
    for i in [0, 1, 2]:
        for j in [0, 1, 2]:
            header.remove(f"A_{i}_{j}", ignore_missing=True)
            header.remove(f"B_{i}_{j}", ignore_missing=True)
            header.remove(f"AP_{i}_{j}", ignore_missing=True)
            header.remove(f"BP_{i}_{j}", ignore_missing=True)
    i = "ORDER"
    header.remove(f"A_{i}", ignore_missing=True)
    header.remove(f"B_{i}", ignore_missing=True)
    header.remove(f"AP_{i}", ignore_missing=True)
    header.remove(f"BP_{i}", ignore_missing=True)
    return header


# fixup headers
def add_NAXES(header, im, add_naxisi=True, inplace=False, verbose=False):
    """Add NAXIS and NAXIS1, NAXIS2 keywords to a header."""
    log("adding NAXES", level=0)
    if not inplace:
        header = header.copy()
    header["NAXIS"] = 2
    if add_naxisi:
        header["NAXIS1"] = im.width
        header["NAXIS2"] = im.height
    return header


def remove_cd(header, verbose=False):
    log("removing CD matrix", level=1)
    header = header.copy()
    if hasattr(WCS(header).wcs, "cd"):
        pre = "CD"
    elif hasattr(WCS(header).wcs, "pc"):
        pre = "PC"
    else:
        return header
    header.remove(pre + "1_1", ignore_missing=True)
    header.remove(pre + "1_2", ignore_missing=True)
    header.remove(pre + "2_1", ignore_missing=True)
    header.remove(pre + "2_2", ignore_missing=True)
    return header


def flip_parity2(header, width=None, inplace=False):
    """
    Flip the parity of the FITS header
    This cannot be done with a header that has only
    scale and rotation matrix, it must have a CD or PC matrix.
    We flip the parity by flipping the 1_2 and 2_2 elements of
    PC/CD matrix
    the height should be from NAXIS2, or given
    """
    if not inplace:
        header = header.copy()
    # original_cd_sign = get_cd_sign(header = header)
    # original_pixel_scale_matrix = get_cd(header)
    if width is None:
        if "NAXIS1" in header:
            width = header["NAXIS1"]
        else:
            raise ValueError("Cannot flip parity without height")

    header["CRPIX2"] = width + 1 - header["CRPIX2"]
    if "PC1_2" in header:
        log("flipping PC", level=0)
        header["PC1_1"] *= -1
        header["PC2_1"] *= -1
    elif "CD1_2" in header:
        log("flipping CD", level=0)
        header["CD1_1"] *= -1
        header["CD2_2"] *= -1
    else:
        raise ValueError("Cannot flip parity without CD or PC matrix")
    # new_cd_sign = get_cd_sign(header = header)
    # new_pixel_scale_matrix = get_cd(header)

    # # pretty print the matrices showing the change
    # log(f'Original CD matrix: cd_sign = {original_cd_sign}')
    # log(pretty_print_matrix(original_pixel_scale_matrix))
    # log(f'New CD matrix: cd_sign = {new_cd_sign}')
    # log(pretty_print_matrix(new_pixel_scale_matrix))

    return header

def blank_header():
    # Create a blank FITS header
    header = Header()
    naxis1 = 256
    naxis2 = 256

    # Set the coordinate system to ICRS
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['CRVAL1'] = 0.0
    header['CRVAL2'] = 0.0
    header['CRPIX1'] = naxis1 / 2 + 0.5
    header['CRPIX2'] = naxis2 / 2 + 0.5
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CD1_1'] = -0.0002777777777777778
    header['CD1_2'] = 0.0
    header['CD2_1'] = 0.0
    header['CD2_2'] = 0.0002777777777777778
    return header