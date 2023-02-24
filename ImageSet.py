from math import sqrt, abs
from astropy.io.fits import Header
from astropy.wcs import WCS

import wcs_helpers as wh

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
    
    def dict(self):
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