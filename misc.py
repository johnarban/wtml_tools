from  aplpy import FITSFigure

from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u 


def preview_image(url,coords=None):
    g =FITSFigure(url)
    g.add_grid()
    g.set_theme('pretty')
    
    if coords is None:
        xy = g._data.shape # get coords at center
        # get coords at center
        ra,dec = g.pixel2world(min(*xy),min(*xy))
        coord_center = SkyCoord(ra=ra,dec=dec,unit=(u.deg, u.deg))
        catalog = Vizier.query_region(coord_center,radius=2*u.deg,catalog='III/135A/catalog') # HD catalog
        coords = catalog[0]['_RA.icrs','_DE.icrs']
        coords.rename_columns(['_RA.icrs','_DE.icrs'],['ra','dec'])
        coords = SkyCoord(coords['ra'],coords['dec'],unit=(u.hourangle, u.deg))

    g.show_markers(coords.ra,coords.dec,marker='x',facecolor='red',edgecolor='red',s=15)
    g.axis_labels.show()
    # g.tick_labels.set_xformat('hh:mm:ss.ss')
    g.tick_labels.show()
    g.tick_labels.set_style('colons')
    g.show_rgb()
    return g