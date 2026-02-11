import os
from importlib import reload
from math import sqrt
import json
from pyavm import AVM
from pyavm.exceptions import NoXMPPacketFound

from astropy.io.fits import Header, ImageHDU
from astropy.wcs import WCS

import io_helpers as ih
import wcs_helpers as wh
import avm_utils as au

import numpy as np

try:
    from anyascii import anyascii
except:
    from string import printable

    anyascii = lambda s: "".join(filter(lambda x: x in set(printable), s))
reload(ih)
reload(wh)
reload(au)

import logger as logger

reload(logger)

FITS_EXTENSIONS = [".fits", ".fit", ".fts", ".fz", ".fits.fz"]
IMAGE_EXTENSIONS = [".jpg", ".jpeg", ".png"]


class ImageHeader:
    def __init__(
        self,
        image_path,
        wcsfile=None,
        image_size=None,
        use_avm=False,
        rotate_cd_matrix=0,
        force_image_center=False,
    ):
        # Initialize attributes
        self.image_path = image_path
        self.basename = os.path.basename(image_path)
        self.dirname = os.path.dirname(image_path)
        self.imagename = os.path.splitext(self.basename)[0]
        self.wcsfile = wcsfile
        self.header = None
        self.credits = ""
        self.description = ""
        self.name = ""
        self.rotate_cd_matrix_by = rotate_cd_matrix
        self.image_size = image_size  # Prioritize this value

        self.use_avm = use_avm
        self.avm = self.get_avm(image_path)

        # Step 1: Use AVM metadata if applicable
        if use_avm and (self.avm is not None):
            self.header = self.header_from_avm(self.avm)
            if self.rotate_cd_matrix_by != 0:
                self.header = wh.rotate_cd_matrix(
                    self.header, self.rotate_cd_matrix_by
                )
            self.check_flip_parity()
        else:
            use_avm = False

            # Step 2: Use WCS file if provided
            if isinstance(self.wcsfile, Header):
                self.header = self.wcsfile
            else:
                self.wcsfile = self._valid_wcs_file(self.wcsfile)
                if self.wcsfile is not None:
                    self.wcs_header = ih.open_header(self.wcsfile)
                    if not self.use_avm:
                        self.header = self.wcs_header

        # Step 3: Create a blank header if none exists
        if self.header is None:
            logger.log(
                f"Could not find header for {self.image_path}", level="ERROR"
            )
            logger.log("Creating blank header", level="ERROR")
            self.header = wh.blank_header()

        # Step 4: Override header dimensions with image_size if provided
        if self.image_size:
            width, height = self.image_size
            self.header["NAXIS1"] = width
            self.header["NAXIS2"] = height
            self.header["IMAGEW"] = width
            self.header["IMAGEH"] = height

        # Step 5: Clean and normalize the header
        self.clean_header()
        self.normalize_header()

        # Step 6: Force CRPIX to center if requested
        if force_image_center:
            self.crpix_to_center()

    def __repr__(self):
        return self.header.__repr__()

    def get_avm(self, filename):
        try:
            if filename.startswith("http"):
                return None
            return AVM.from_image(filename)
        except NoXMPPacketFound:
            return None

    def _valid_wcs_file(self, wcsfile):
        """
        Check if the WCS file exists and return the path to the file.
        If the WCS file is not found, try to find it in the image directory.
        If the WCS file is not found, return None.
        """
        if wcsfile is None:
            wcsfile = ih.guess_wcs_filename(self.image_path)
            if wcsfile is None:
                logger.log(
                    f"Could not find WCS file for {self.image_path}",
                    level="ERROR",
                )
                raise Exception(
                    f"Could not find WCS file for {self.image_path}"
                )
            else:
                return wcsfile
        else:
            wcs_dirname = os.path.dirname(wcsfile)
            wcs_basename = os.path.basename(wcsfile)

            if wcs_basename.find("*") == 0:
                wcs_basename = wcs_basename.replace(
                    "*", os.path.splitext(self.imagename)[0]
                )
                wcsfile = os.path.join(wcs_dirname, wcs_basename)

            if os.path.exists(wcsfile):
                return wcsfile

            elif os.path.exists(os.path.join(self.dirname, wcs_basename)):
                logger.log(
                    f"Found WCS file in image directory (not on {wcs_dirname})",
                    "INFO",
                )
                return os.path.join(self.dirname, wcs_basename + ".wcs")
            else:
                logger.log(
                    f"Could not find WCS file {wcsfile} for {self.image_path}",
                    level="ERROR",
                )
                raise Exception(
                    f"Could not find WCS file {wcsfile} for {self.image_path}"
                )

    def clean_header(self):
        """
        Remove all keys from the header that are not needed for the WTML file.
        """
        self.header = wh.clean_header(self.header)
        self.header = wh.remove_sip(self.header)

    def normalize_header(self):
        """
        Normalize the header to ensure consistency in dimensions and required keywords.
        """
        header = self.header.copy()
        header = wh.add_cd_matrix(header)

        # Ensure NAXIS1, NAXIS2, IMAGEW, and IMAGEH are consistent
        par = self._get_naxis_params()
        header = wh.add_NAXES(header, **par)

        # Add required FITS keywords
        header = wh.add_required_keywords(header)

        self.header = header

    def crpix_to_center(self):
        """
        Move CRPIX to the center of the image based on the dimensions.
        """
        logger.log("Moving CRPIX to center", level="INFO")
        width, height = self._get_naxis_params().values()
        header = self.header.copy()
        header["CRPIX1"] = (width + 1) / 2
        header["CRPIX2"] = (height + 1) / 2
        self.header = header

    def check_flip_parity(self):
        par = self._get_naxis_params()

        header = wh.add_cd_matrix(self.header)

        # check parity
        if not wh.is_JPEGLike(parity=wh.get_parity(header=header)):
            logger.log(
                "WTML requires JPEG-like parity (parity < 0, positive det(CD))",
                level="ERROR",
            )
            logger.log("Flipping parity", level="ERROR")
            self.header = wh.flip_parity(header, par["height"])

    def _get_naxis_params(self):
        """
        Get the image dimensions (width and height) from image_size, header, or file.
        """
        if self.image_size is not None:
            width, height = self.image_size
        elif self.header is not None and "IMAGEH" in self.header:
            width, height = self.header["IMAGEW"], self.header["IMAGEH"]
        else:
            width, height = ih.get_image_size(self.image_path)
        return {"width": width, "height": height}

    def add_naxis_params(self):
        """
        Add NAXIS1 and NAXIS2 to the header.
        """
        par = self._get_naxis_params()
        self.header = wh.add_NAXES(header=self.header, **par)

    @property
    def rotated_ccd_header(self):
        return self.header

    def header_from_avm(self, avm):
        par = self._get_naxis_params()
        if not hasattr(avm, "Spatial"):
            return None

        header = avm.to_wcs().to_header()

        header["IMAGEW"] = par["width"]
        # add comment
        header.comments["IMAGEW"] = "Image width in pixels"
        header["IMAGEH"] = par["height"]
        header.comments["IMAGEH"] = "Image height in pixels"

        if "Description" in avm._items:
            header["DESCRIPT"] = anyascii(avm._items["Description"] or "none")
            header.comments["DESCRIPT"] = "Description of the image"
        if "Credits" in avm._items:
            header["CREDITS"] = avm._items["Credits"]
            header.comments["CREDITS"] = "Credits for the image"

        return header

    def write_avm(self, path=".", force_180=False):
        image_name = os.path.basename(self.image_path)
        au.write_avm(
            image_name,
            self.header,
            suffix="_avm",
            path_out=path,
            force_180=force_180,
        )

    def write_fits(self, path=".", name=None):
        data = ih.get_PIL_image(self.image_path)
        # convert PIL image to grayscale
        if "RGB" in data.mode:
            data = data.convert("L")
        if name is None:
            name = self.imagename + ".fits"
        fits_path = os.path.join(path, name)
        # output image as fits files

        hdu = ImageHDU(data=data, header=self.header)
        hdu.writeto(fits_path, overwrite=True)

    def write_header(self, name=None):
        if name is None:
            name = self.imagename + ".wcs.fits"
        header_path = os.path.join(self.dirname, name)
        # output header as fits files
        header = wh.fixup_header(self.header)
        header.tofile(
            header_path, sep="", endcard=True, padding=True, overwrite=True
        )


class ImageSet:

    def __init__(
        self,
        DemUrl="",
        MSRCommunityId="0",
        MSRComponentId="0",
        Permission="0",
        Generic="False",
        DataSetType="Sky",
        BandPass="Visible",
        Url="",
        TileLevels="0",
        WidthFactor="2",
        Sparse="True",
        Rotation="0",
        QuadTreeMap="",
        Projection="SkyImage",
        Name="",
        FileType=".jpg",
        CenterY=0,
        CenterX=0,
        BottomsUp=None,
        parity=None,
        StockSet=False,
        ElevationModel=False,
        OffsetX=0,
        OffsetY=0,
        BaseTileLevel=0,
        BaseDegreesPerTile=1,
        ReferenceFrame="",
        MeanRadius="0",
    ):

        self.DemUrl = DemUrl
        self.MSRCommunityId = MSRCommunityId
        self.MSRComponentId = MSRComponentId
        self.StockSet = StockSet
        self.Permission = Permission
        self.Generic = Generic

        self.DataSetType = DataSetType
        self.ElevationModel = ElevationModel
        self.QuadTreeMap = QuadTreeMap
        self.Projection = Projection
        self.Sparse = Sparse

        self.BandPass = BandPass

        self.Name = Name

        self.Url = Url
        self.FileType = FileType

        self.TileLevels = TileLevels
        self.WidthFactor = WidthFactor
        self.BaseTileLevel = BaseTileLevel
        self.BaseDegreesPerTile = BaseDegreesPerTile
        self.MeanRadius = MeanRadius

        self.Rotation = Rotation

        self.BottomsUp = BottomsUp
        if self.BottomsUp is None:
            if parity is not None:
                self.BottomsUp = parity > 0
            else:
                self.BottomsUp = False
        self.CenterY = CenterY
        self.CenterX = CenterX
        self.OffsetX = OffsetX
        self.OffsetY = OffsetY

        self.ReferenceFrame = ReferenceFrame

    def set_wtml_params_from_header(self, header):
        header = header.copy()

        wcs = WCS(header)
        crpix = [(header["NAXIS1"] + 1) / 2, (header["NAXIS2"] + 1) / 2]
        crval = wcs.wcs_pix2world(*crpix, 1)
        offset = crpix
        scales, rot, parity = wh.get_scale_rot(header)
        scale = sqrt(abs(scales[0] * scales[1]))
        bottoms_up = parity > 0

        self.CenterX = crval[0]
        self.CenterY = crval[1]
        self.OffsetX = offset[0]
        self.OffsetY = offset[1]
        self.BottomsUp = bottoms_up
        self.Rotation = rot
        if self.Projection == "SkyImage":
            self.BaseDegreesPerTile = scale
        else:
            raise Exception("Projection not supported")

    def json(self):
        return {k: str(v) for k, v in self.__dict__.items()}

    def set_name(self, name):
        self.Name = name

    def set_url(self, url):
        self.Url = url

    def set_thumbnail(self, thumbnail):
        self.ThumbnailURL = thumbnail

    def set_projection(self, projection):
        self.Projection = projection

    def set_tile_levels(self, levels):
        self.TileLevels = levels

    def set_base_tile_level(self, level):
        self.BaseTileLevel = level

    def set_base_degrees_per_tile(self, degrees):
        self.BaseDegreesPerTile = degrees

    def set_bottoms_up(self, bottoms_up):
        self.BottomsUp = bottoms_up
