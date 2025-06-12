from astropy.io.fits import Header, PrimaryHDU
from astropy.wcs import WCS
from math import atan2, atan
from astropy.wcs.utils import _is_cd_orthogonal
from numpy import pi
import logger as logger
from importlib import reload

reload(logger)
import numpy as np


def ensure_wcs(supposed_wcs):
    """Converts a header to a wcs object if it isn't already one."""
    if isinstance(supposed_wcs, WCS):
        return supposed_wcs
    else:
        return WCS(supposed_wcs)


def add_required_keywords(header, inplace=False):
    # need to ensure SIMPLE and BITPIX are in the header
    if not inplace:
        header = header.copy()
    if "SIMPLE" not in header:
        # add simple as first card
        header.insert(0, ("SIMPLE", True))
        # header['SIMPLE'] = True
    if "BITPIX" not in header:
        # add bitpix as second card
        header.insert(1, ("BITPIX", -32))
        # header['BITPIX'] = -32
    # and end card to make sure it is a valid FITS file
    return header


def fixup_header(header):
    # convert a header to a wcs and back to a header
    # and add back in any missing keywords
    # this is useful for cleaning up a header
    wcs = ensure_wcs(header)
    new_header = wcs.to_header()
    new_header = PrimaryHDU(header=new_header).header
    # new_header = add_required_keywords(new_header)
    # # # fill in keys that are in header but missing in new_header
    # for key in header:
    #     if (
    #         (key not in new_header)
    #         and (key not in ['SIMPLE', 'BITPIX'])
    #         and not ('CD' in key and '_' in key)
    #         and not ('PC' in key and '_' in key)
    #         ):
    #         print(key)
    #         new_header[key] = header[key]
    # add back naxis1 and naxis2, naxis if present
    if "NAXIS" in header:
        new_header["NAXIS"] = header["NAXIS"]
    if "NAXIS1" in header:
        new_header["NAXIS1"] = header["NAXIS1"]
    if "NAXIS2" in header:
        new_header["NAXIS2"] = header["NAXIS2"]
    # imagew and imageh
    if "IMAGEW" in header:
        new_header["IMAGEW"] = header["IMAGEW"]
    if "IMAGEH" in header:
        new_header["IMAGEH"] = header["IMAGEH"]

    # replace pc matrix with cd matrix and cdelt
    # new_header = convert_pc_to_cd(new_header)

    return new_header


def convert_pc_to_cd(header, inplace=False):
    cd_matrix = get_cd(header)
    wcs = ensure_wcs(header)
    cdelt = wcs.wcs.get_cdelt()
    if not inplace:
        header = header.copy()
    header["CD1_1"] = cd_matrix[0, 0]
    header["CD1_2"] = cd_matrix[0, 1]
    header["CD2_1"] = cd_matrix[1, 0]
    header["CD2_2"] = cd_matrix[1, 1]
    header["CDELT1"] = cdelt[0]
    header["CDELT2"] = cdelt[1]
    header.remove("PC1_1", ignore_missing=True)
    header.remove("PC1_2", ignore_missing=True)
    header.remove("PC2_1", ignore_missing=True)
    header.remove("PC2_2", ignore_missing=True)

    return header


def ensure_header(supposed_header):
    """Converts a wcs to a header if it isn't already one."""
    if isinstance(supposed_header, Header):
        return supposed_header
    else:
        return supposed_header.to_header()


def is_header_or_wcs(supposed_wcs_or_header):
    """Returns True if the object is a header or a wcs object."""
    return isinstance(supposed_wcs_or_header, Header) or isinstance(
        supposed_wcs_or_header, WCS
    )


## CD Matrix tools
## FITS originally used CDELT and CROTA to describe the scale and rotation
## Then it switched to using
def get_cd(header=None, wcs=None, as_dict=False):
    """Return the CD matrix from a header or wcs object.
    This is a convenience function that will return the CD matrix
    from either a header or a wcs object.
    This works using the astropy.wcs.WCS.pixel_scale_matrix property.
    This property derives from either wcs.wcs.cd or wcs.wcs.pc & wcs.wcs.cdelt.
    """
    # check if user was lazy and header is actually a wcs object
    if header is not None:
        wcs = ensure_wcs(header)
    elif wcs is not None:
        wcs = ensure_wcs(wcs)

    if as_dict:
        cd = wcs.pixel_scale_matrix
        return {
            "CD1_1": cd[0, 0],
            "CD1_2": cd[0, 1],
            "CD2_1": cd[1, 0],
            "CD2_2": cd[1, 1],
        }

    return wcs.pixel_scale_matrix


def det_cd(cd_matrix):
    """Return the determinant of a CD matrix.

    the cd Matrix looks like [ CD1_1 CD1_2 ]
                             [ CD2_1 CD2_2 ]
    written in the form of a numpy array:
                             [ [0,0] [0,1] ]
                             [ [1,0] [1,1] ]
    """
    return cd_matrix[0, 0] * cd_matrix[1, 1] - cd_matrix[0, 1] * cd_matrix[1, 0]


def add_cd_matrix(header):
    """
    add cd matrix to header if one does not exist
    """
    cd_matrix = get_cd(header=header)
    # check if keys starting with CD exist
    matrix_keys = {
        "CD1_1": (0, 0),
        "CD1_2": (0, 1),
        "CD2_1": (1, 0),
        "CD2_2": (1, 1),
    }
    if any([key in header for key in matrix_keys]):
        return header

    for key, value in matrix_keys.items():
        header[key] = cd_matrix[value]
    # remove CDELT and CROTA
    # header.remove('CDELT1', ignore_missing=True)
    # header.remove('CDELT2', ignore_missing=True)
    # header.remove('CROTA2', ignore_missing=True)

    return header


def replace_cd_matrix(header, cd_matrix):
    """
    replace cd matrix in header with cd_matrix
    """
    # check if keys starting with CD exist
    matrix_keys = {
        "CD1_1": (0, 0),
        "CD1_2": (0, 1),
        "CD2_1": (1, 0),
        "CD2_2": (1, 1),
    }
    for key, value in matrix_keys.items():
        header[key] = cd_matrix[value]

    return header


def replace_pc_matrix(header, pc_matrix):
    """
    replace pc matrix in header with pc_matrix
    """
    # check if keys starting with CD exist
    matrix_keys = {
        "PC1_1": (0, 0),
        "PC1_2": (0, 1),
        "PC2_1": (1, 0),
        "PC2_2": (1, 1),
    }
    for key, value in matrix_keys.items():
        header[key] = pc_matrix[value]

    return header


def pretty_print_matrix(matrix, name=""):
    """Prints a matrix in a pretty way."""
    str1 = f" [ {matrix[0,0]:9.3g} {matrix[0,1]:9.3g} ]"
    str2 = f" [ {matrix[1,0]:9.3g} {matrix[1,1]:9.3g} ]"
    return str1 + "\n" + str2


def get_cd_sign(cd=None, header=None):
    """get's the cd_sign of the cd matrix
    FITS typically have positive cd_sign
    while JPEGS tend to have negative cd_sign
    See https://github.com/WorldWideTelescope/toasty/blob/master/toasty/builder.py#L290
    for a discussion on cd_sign
    """

    if (cd is None) and (header is None):
        raise ValueError("Must pass either a cd matrix or a header")

    if cd is None:
        if is_header_or_wcs(header):
            cd = get_cd(wcs=ensure_wcs(header))
    else:
        if is_header_or_wcs(cd):
            logger.log(f"Passed {type(cd)} as cd matrix", level="DEBUG")
            cd = get_cd(header=cd)
        elif is_header_or_wcs(header):
            cd = get_cd(header=header)
    # else cd must be a numpy array

    if det_cd(cd) >= 0:
        return 1
    else:
        return -1


def center_header(header, loc=None, inplace=False):
    if loc is None:
        wcs = WCS(header)
        crpix = [(header["NAXIS1"] + 1) / 2, (header["NAXIS2"] + 1) / 2]
        crval = np.atleast_1d(wcs.wcs_pix2world(*crpix, 1))
    else:
        crpix = loc["crpix"]
        crval = loc["crval"]

    if not inplace:
        header = header.copy()

    header["CRPIX1"] = crpix[0]
    header["CRPIX2"] = crpix[1]
    header["CRVAL1"] = crval[0]
    header["CRVAL2"] = crval[1]

    return header


def modify_header(
    header,
    scale=0,
    rot=0,
    nudge={"x": 0, "y": 0},
    center=None,
    swap=False,
    flipx=False,
    flipy=False,
    inplace=False,
    verbose=False,
):

    if not inplace:
        header = header.copy()
    # scale by modifying CDELT
    header = fixup_header(
        header
    )  # converts to regular form with PC matrix and CDELT = 1
    if scale != 0:
        logger.log(f"Scaling by {scale}", level="ALWAYS")
        if isinstance(scale, dict):
            header["CDELT1"] *= scale["x"]
            header["CDELT2"] *= scale["y"]
        else:
            header["CDELT1"] *= scale
            header["CDELT2"] *= scale

    # rotate the pc matrix
    if rot != 0:
        logger.log(f"Rotating by {rot} degrees", level="ALWAYS")
        header = rotate_pc_matrix(header, angle=rot, ccw=True, swap=swap)

    # makes sure the rotation is about the center of the imag

    if center is not None:
        if "NAXIS1" in header and "NAXIS2" in header:
            if isinstance(center, dict):
                header = center_header(header, loc=center)
            else:
                header = center_header(header)
        else:
            logger.log(
                "Can't center image if nor naxis keywords are present",
                level="ERROR",
            )

    if nudge["x"] != 0 or nudge["y"] != 0:
        logger.log(f"Nudging by {nudge}", level="ALWAYS")
        # nudge the crpix
        header["CRPIX1"] += nudge["x"]
        header["CRPIX2"] += nudge["y"]

    if flipy:
        logger.log("Flipping Y axis", level="ALWAYS")
        header["CDELT2"] *= -1
    if flipx:
        logger.log("Flipping X axis", level="ALWAYS")
        header["CDELT1"] *= -1

    return header


def get_parity(cd=None, header=None):
    """Return the parity of the CD matrix."""
    return -get_cd_sign(cd=cd, header=header)


def is_JPEGLike(parity=None, cd_sign=None):
    """Return True if the parity is negative."""
    if is_header_or_wcs(parity):
        logger.log(
            "Passed a header or wcs object to is_JPEGLike, assuming it is a header",
            level="DEBUG",
        )
        parity = get_parity(header=parity)
    if cd_sign is not None:
        return cd_sign > 0
    return parity < 0


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
            raise Exception("Cannot flip parity without height")

    header["CRPIX2"] = height + 1 - header["CRPIX2"]
    if "PC1_2" in header:
        logger.log("flipping PC", level="DEBUG")
        header["PC1_2"] *= -1
        header["PC2_2"] *= -1
    elif "CD1_2" in header:
        logger.log("flipping CD", level="DEBUG")
        header["CD1_2"] *= -1
        header["CD2_2"] *= -1
    else:
        raise Exception("Cannot flip parity without CD or PC matrix")
    new_cd_sign = get_cd_sign(header=header)
    new_pixel_scale_matrix = get_cd(header)

    # pretty print the matrices showing the change
    logger.log(
        f"Original CD matrix: cd_sign = {original_cd_sign}", level="INFO"
    )
    logger.log(pretty_print_matrix(original_pixel_scale_matrix), level="DEBUG")
    logger.log(f"New CD matrix: cd_sign = {new_cd_sign}", level="INFO")
    logger.log(pretty_print_matrix(new_pixel_scale_matrix), level="DEBUG")

    return header


def flip_parity2(header, width=None, inplace=False):
    """
    Flip the parity of the FITS header (north-south flip)
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
    if width is None:
        if "NAXIS1" in header:
            width = header["NAXIS1"]
        else:
            raise ValueError("Cannot flip parity without width")

    header["CRPIX1"] = width + 1 - header["CRPIX1"]
    if "PC1_2" in header:
        logger.log("flipping PC", level="DEBUG")
        header["PC1_1"] *= -1
        header["PC2_1"] *= -1
    elif "CD1_2" in header:
        logger.log("flipping CD", level="DEBUG")
        header["CD1_1"] *= -1
        header["CD2_1"] *= -1
    else:
        raise ValueError("Cannot flip parity without CD or PC matrix")
    new_cd_sign = get_cd_sign(header=header)
    new_pixel_scale_matrix = get_cd(header)

    # pretty print the matrices showing the change
    logger.log(
        f"Original CD matrix: cd_sign = {original_cd_sign}", level="INFO"
    )
    logger.log(pretty_print_matrix(original_pixel_scale_matrix), level="DEBUG")
    logger.log(f"New CD matrix: cd_sign = {new_cd_sign}", level="INFO")
    logger.log(pretty_print_matrix(new_pixel_scale_matrix), level="DEBUG")

    return header


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
    logger.log("Finding the rotation angle", level="DEBUG")
    cd_sign = get_cd_sign(cd)
    logger.log(f"\tCD sign ={cd_sign}", level="DEBUG")
    # wwt
    T = cd[1, 1]
    A = cd[0, 1]
    rot_wwt = atan2(-cd_sign * A, -T) * 180 / pi
    logger.log(f"\tRot (WWT) = {rot_wwt:0.3f} degrees", level="DEBUG")

    # astrometry.net
    T = cd_sign * cd[0, 0] + cd[1, 1]
    A = cd_sign * cd[1, 0] - cd[0, 1]
    rot_astrom = -atan2(-cd_sign * A, -T) * 180 / pi
    logger.log(
        f"\tRot (Astrometry.net) = {rot_astrom:0.3f} degrees", level="DEBUG"
    )
    if wwt:
        logger.log(f"\tUsing WWT method. Angle :{rot_wwt:0.3f}", level="DEBUG")
        return rot_wwt
    else:
        logger.log(
            f"\tUsing Astrometry.net method. Angle :{rot_astrom:0.3f}",
            level="DEBUG",
        )
        return rot_astrom


def get_scale(cd):
    """Return the scale of the CD matrix."""
    return (cd**2).sum(axis=0) ** 0.5


def get_scale_rot(header, wwt=False):
    # from astrometry.net/net/wcs.py
    wcs = ensure_wcs(header)
    cd = get_cd(wcs=wcs)

    scales = get_scale(cd)  # in degrees / pixel

    rot = get_rot_from_cd(cd, wwt=wwt)  # in degress

    parity = get_parity(cd)  # 1 or -1

    logger.log("get_scale_rot::", level="DEBUG")
    logger.log(
        f"\tscale: {scales}\n\trot: {rot:.3f}\n\tparity: {parity}",
        level="DEBUG",
    )

    return scales, rot, parity


def clean_header(header, inplace=False, verbose=False):
    """Remove all HISTORY and COMMENT cards from a header."""
    logger.log("cleaning header", level="DEBUG")
    if not inplace:
        header = header.copy()
    header.remove("HISTORY", remove_all=True, ignore_missing=True)
    header.remove("COMMENT", remove_all=True, ignore_missing=True)
    return header


def rotation_matrix(angle, radians=False):
    angle = np.radians(angle) if not radians else angle
    return np.array(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )


def rotate_cd_matrix(header, angle=90, ccw=True, swap=True):
    if "CD1_1" not in header:
        if "PC1_1" in header:
            logger.log(
                "CD matrix not present. Attempting to rotate PC matrix",
                level="INFO",
            )
            header = rotate_pc_matrix(header, angle=angle, ccw=ccw)
        else:
            logger.log(
                "Cannot rotate CD matrix without CD or PC matrix", "ERROR"
            )
            return header
    # rotate CD matrix by 90 degrees (counter-clock-wise)
    logger.log(f"Rotating CD Matrix by {angle} degrees", level="INFO")

    hdr = header.copy()

    cd_matrix = get_cd(header)
    if ccw:
        cd_matrix = rotation_matrix(angle) @ cd_matrix
    else:
        cd_matrix = rotation_matrix(-angle) @ cd_matrix

    hdr = replace_cd_matrix(hdr, cd_matrix)

    if swap:
        # swap CRPIX1 and CRPIX2
        logger.log("\talso swapping CRPIX1 and CRPIX2", level="INFO")
        hdr["CRPIX1"], hdr["CRPIX2"] = header["CRPIX2"], header["CRPIX1"]

    return hdr


def get_pc(header):
    wcs = WCS(ensure_header(header))
    return wcs.wcs.get_pc()


def rotate_pc_matrix(header, angle=90, ccw=True, swap=True):
    if "PC1_1" not in header:
        if "CD1_1" in header:
            logger.log(
                "PC matrix not present. Attempting to rotate CD matrix",
                level="INFO",
            )
            header = rotate_cd_matrix(header, angle=angle, ccw=ccw)
        else:
            logger.log(
                "Cannot rotate PC matrix without CD or PC matrix", "ERROR"
            )
            return header

    logger.log(f"Rotating PC Matrix by {angle} degrees", level="INFO")

    hdr = header.copy()

    pc_matrix = get_pc(header)
    if ccw:
        pc_matrix = rotation_matrix(angle) @ pc_matrix
    else:
        pc_matrix = rotation_matrix(-angle) @ pc_matrix

    hdr = replace_pc_matrix(hdr, pc_matrix)

    if swap:
        # swap CRPIX1 and CRPIX2
        logger.log("\talso swapping CRPIX1 and CRPIX2", level="INFO")
        hdr["CRPIX1"], hdr["CRPIX2"] = header["CRPIX2"], header["CRPIX1"]

    return hdr


def remove_sip(header, inplace=False, verbose=False):
    """Remove all SIP related keywords from a header."""
    if verbose:
        logger.log("removing SIP", level="DEBUG")
    if not inplace:
        header = header.copy()
    if "CTYPE1" in header:
        header["CTYPE1"] = header["CTYPE1"].replace("-SIP", "")
    if "CTYPE2" in header:
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
def add_NAXES(header, width=None, height=None, inplace=False, verbose=False):
    """Add NAXIS and NAXIS1, NAXIS2 keywords to a header."""
    logger.log("adding NAXES", level="DEBUG")
    if not inplace:
        header = header.copy()
    header["NAXIS"] = 2
    if width is not None:
        header["NAXIS1"] = width
    elif "IMAGEW" in header:
        header["NAXIS1"] = header["IMAGEW"]
    if height is not None:
        header["NAXIS2"] = height
    elif "IMAGEH" in header:
        header["NAXIS2"] = header["IMAGEH"]

    return header


def remove_cd(header, verbose=False):
    logger.log("removing CD matrix", level="DEBUG")
    header = header.copy()

    pre = "CD"
    header.remove(pre + "1_1", ignore_missing=True)
    header.remove(pre + "1_2", ignore_missing=True)
    header.remove(pre + "2_1", ignore_missing=True)
    header.remove(pre + "2_2", ignore_missing=True)

    pre = "PC"
    header.remove(pre + "1_1", ignore_missing=True)
    header.remove(pre + "1_2", ignore_missing=True)
    header.remove(pre + "2_1", ignore_missing=True)
    header.remove(pre + "2_2", ignore_missing=True)
    return header


def blank_header():
    # Create a blank FITS header
    header = Header({"SIMPLE": True, "BITPIX": -32})
    naxis1 = 256
    naxis2 = 256

    # Set the coordinate system to ICRS
    header["BITPIX"] = -32
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    header["CRVAL1"] = 0.0
    header["CRVAL2"] = 0.0
    header["CRPIX1"] = naxis1 / 2 + 0.5
    header["CRPIX2"] = naxis2 / 2 + 0.5
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["CD1_1"] = -10 / 3600
    header["CD1_2"] = 0.0
    header["CD2_1"] = 0.0
    header["CD2_2"] = 10 / 3600
    return header
