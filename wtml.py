from io_helpers import guess_wcs_filename
class WTML:
    
    def __init__(self, image_path, wcs_path = None):
        self.header = None
        self.wcs = None
        self.image_path = None
        self.wcs_path = None
    
    def from_image(self, image_path, wcs_path):
        return None
    
    def from_avm(self, image_path):
        return None
    
    def from_fits(self, image_path):
        return None
        