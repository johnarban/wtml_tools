    
# no currently implemented
class WTMLHeader:
    
    def __init__(self, image_path, wcs_file=None):
        
        self.image_path = image_path
        self.wcs_file = wcs_file
        
        self.wcs = None
        self.header = None
        self.wtml = None
        
        self._load_wcs()
        self._load_header()
        self._load_wtml()
        
    
    @classmethod
    def from_avm(cls, avm):
        avm_header = avm.to_wcs().to_header()
        avm_header["NAXIS"] = 2
        avm_header["NAXIS1"] = avm.Spatial.ReferenceDimension[0]
        avm_header["NAXIS2"] = avm.Spatial.ReferenceDimension[1]
        avm_header = wh.flip_parity(avm_header)
        return cls.from_header(avm_header)
    
    @classmethod
    def from_header(cls, header):
        return cls.from_wcs(WCS(header))
    
    @classmethod
    def from_wcs(cls, wcs):
        return cls.from_wcs(wcs)

from wcs_helpers import get_scale_rot, get_parity, is_JPEGLike, flip_parity, remove_cd, blank_header

class TaggedImage:
    
    def __init__(self, image_path, wcs_file=None):
        
        self.image_path = image_path
        
        self.wcs_file = wcs_file
        
        self.filename = os.path.basename(image_path)
        self.directory = os.path.dirname(image_path)
        
        self.wwt_url = None
        
        self._has_avm = False # if image is already avm tagged we don't need to write out a tagged image
        self.avm = None
        
        self._get_image_header()
        
        if self.avm is None:
            log("No AVM deteected. Creating AVM from self.header", level='INFO')
            self.avm = AVM.from_header(self.header)
        
        
    def _get_image_header(self):
        """
        Get the image header from a FITS file.
        """
        if self.wcs_file is not None:
            return Header.fromfile(self.wcs_file)
        else:
            if os.path.splitext(self.image_file)[1] in FITS_EXTENSIONS:
                log("Reading FITS header from image file", level=1)
                return Header.fromfile(self.image_file), None
            elif os.path.splitext(self.image_file)[1] in IMAGE_EXTENSIONS:
                try:
                    log("Reading AVM from image file", level=1)
                    self.avm = AVM.from_file(self.image_file)
                    self._has_avm = True
                    self.header = self.avm.to_wcs().to_header()
                except:
                    log("Could not read AVM from image file", level=2)
                    log("Returning blank header")
                    return blank_header(), None
    
    
    def _get_wtml(self):
        
        
    