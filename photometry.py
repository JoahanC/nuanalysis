import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from shapely.plotting import plot_polygon
from photutils.aperture import RectangularAperture, CircularAperture, CircularAnnulus
import matplotlib
from astropy import units as u



def calculate_background_area(data, src_x, src_y, optimal_radius, display=False):
    """
    Calculates the background area of a circular annulus defined within a 
    NuSTAR fits image within DET1 coordinates

    Arguments:
    ----------
    data : 
    src_x : float
    src_y : float
    optimal_radius : float
    display : bool
    

    Returns:
    --------

    """

    inner_circle = Point(src_x, src_y).buffer(optimal_radius, resolution=1000)
    outer_circle = Point(src_x, src_y).buffer(optimal_radius + 20, resolution=1000)
    horizontal_coords = ((15., 180.), (15., 182.), (345., 182.), (345., 180.)) 
    inner_horizontal_bar = Polygon(horizontal_coords)
    vertical_coords = ((180., 15.), (182., 15.), (182., 345.), (180., 345.)) 
    inner_vertical_bar = Polygon(vertical_coords)
    difference = outer_circle.difference(inner_circle)
    frame_coords = ((15., 15.), (15., 345.), (345., 345.), (345., 15.))
    rectangle = Polygon(frame_coords)
    shared_region = difference.intersection(rectangle)
    final_region = shared_region.difference(inner_horizontal_bar).difference(inner_vertical_bar)

    #fig, ax = plt.subplots()
    #ax.imshow(data, norm=matplotlib.colors.LogNorm())
    #plot_polygon(final_region, ax=ax, add_points=False)
    #plt.show()
    return final_region.area


def calculate_source_area(data, src_x, src_y, optimal_radius):

    circle = Point(src_x, src_y).buffer(optimal_radius, resolution=1000)
    horizontal_coords = ((15., 180.), (15., 182.), (345., 182.), (345., 180.)) 
    inner_horizontal_bar = Polygon(horizontal_coords)
    vertical_coords = ((180., 15.), (182., 15.), (182., 345.), (180., 345.)) 
    inner_vertical_bar = Polygon(vertical_coords)
    frame_coords = ((15., 15.), (15., 345.), (345., 345.), (345., 15.))
    rectangle = Polygon(frame_coords)
    shared_region = circle.intersection(rectangle)
    final_region = shared_region.difference(inner_horizontal_bar).difference(inner_vertical_bar)
    aper = CircularAperture([src_x, src_y], optimal_radius)
    src_counts, cts_err = aper.do_photometry(data)

    #fig, ax = plt.subplots()
    #ax.imshow(data, norm=matplotlib.colors.LogNorm())
    #plot_polygon(final_region, ax=ax, add_points=False)
    #aper.plot(ax)
    #plt.show()
    return final_region.area

def source_counts(data, src_x, src_y, optimal_radius):

    aper = CircularAperture([src_x, src_y], optimal_radius)
    src_counts, cts_err = aper.do_photometry(data)
    #fig, ax = plt.subplots()
    #ax.imshow(data, norm=matplotlib.colors.LogNorm())
    #aper.plot(ax)
    #plt.show()
    #plt.close()
    return src_counts[0] * u.ct

def bkg_counts(data, src_x, src_y, optimal_radius):

    aper = CircularAnnulus([src_x, src_y], optimal_radius, optimal_radius + 20)
    bkg_counts, cts_err = aper.do_photometry(data)
    #fig, ax = plt.subplots()
    #ax.imshow(data, norm=matplotlib.colors.LogNorm())
    #aper.plot(ax)
    #plt.show()
    #plt.close()
    return bkg_counts[0] * u.ct