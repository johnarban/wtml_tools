import os
from importlib import reload
from math import sqrt
import json

from astropy.io.fits import Header
from astropy.wcs import WCS

import io_helpers as ih
import wcs_helpers as wh

reload(ih)
reload(wh)

DEBUG_LEVELS = {"DEBUG": 0, "INFO": 1}
DEBUG_LEVEL = 1
FITS_EXTENSIONS = [".fits", ".fit", ".fts", ".fz", ".fits.fz"]
IMAGE_EXTENSIONS = [".jpg", ".jpeg", ".png"]


def set_debug_level(level):
    global DEBUG_LEVEL
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


class ImageHeader:
    def __init__(self, image_path, wcsfile=None, use_avm=False):
        self.image_path = image_path
        self.basename = os.path.basename(image_path)
        self.dirname = os.path.dirname(image_path)
        self.imagename = os.path.splitext(self.basename)[0]

        self.header = None

        self.use_avm = use_avm
        self.avm = ih.get_avm(image_path)
        if use_avm and self.avm is not None:
            self.header = self.avm.to_wcs().to_header()
        else:
            use_avm = False

        self.wcsfile = self._valid_wcs_file(wcsfile)
        if self.wcsfile is not None:
            self.wcs_header = Header.fromfile(self.wcsfile)
            if not self.use_avm:
                self.header = self.wcs_header

        if self.header is None:
            log(f"Could not find header for {self.image_path}", level="INFO")
            log(f"Creating blank header", level="INFO")
            self.header = wh.blank_header()

    def _valid_wcs_file(self, wcsfile):
        """
        Check if the WCS file exists and return the path to the file.
        If the WCS file is not found, try to find it in the image directory.
        If the WCS file is not found, return None.
        """
        if wcsfile is None:
            wcsfile = ih.guess_wcs_filename(self.image_path)
            if wcsfile is None:
                log(f"Could not find WCS file for {self.image_path}", level="INFO")
                raise Exception(f"Could not find WCS file for {self.image_path}")
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
                log(f"Found WCS file in image directory (not on {wcs_dirname})", "INFO")
                return os.path.join(self.dirname, wcs_basename + ".wcs")
            else:
                log(
                    f"Could not find WCS file {wcsfile} for {self.image_path}",
                    level="INFO",
                )
                raise Exception(
                    f"Could not find WCS file {wcsfile} for {self.image_path}"
                )




class ImageSet:
    
    
    def __init__(self, 
                DemUrl = "",
                MSRCommunityId = "0",
                MSRComponentId = "0",
                Permission = "0",
                Generic = "False",
                DataSetType = "Sky",
                BandPass = "Visible",
                Url = "",
                TileLevels = "0",
                WidthFactor = "2",
                Sparse = "True",
                Rotation = "0",
                QuadTreeMap = "",
                Projection = "SkyImage",
                Name = "",
                FileType = ".jpg",
                CenterY = 0,
                CenterX = 0,
                BottomsUp = False,
                StockSet = False,
                ElevationModel = False, 
                OffsetX = 0,
                OffsetY = 0,
                BaseTileLevel = 0,
                BaseDegreesPerTile = 1,
                ReferenceFrame = "",
                MeanRadius = "0"
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
        return {k:str(v) for k,v in self.__dict__.items()}
        
    def set_name(self,name):
        self.Name = name
        
    def set_url(self,url):
        self.Url = url
    
    def set_thumbnail(self,thumbnail):
        self.ThumbnailURL = thumbnail
    
    def set_projection(self,projection):
        self.Projection = projection
    
    def set_tile_levels(self,levels):
        self.TileLevels = levels
    
    def set_base_tile_level(self,level):
        self.BaseTileLevel = level
    
    def set_base_degrees_per_tile(self,degrees):
        self.BaseDegreesPerTile = degrees
    
    def set_bottoms_up(self,bottoms_up):
        self.BottomsUp = bottoms_up