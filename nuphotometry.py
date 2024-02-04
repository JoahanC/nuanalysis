"""
Boilerplate code for the background estimate routine used in nuanalysis.
"""
import matplotlib
import matplotlib.pyplot as plt
from shapely.plotting import plot_polygon
from shapely.geometry import Point, Polygon
from photutils.aperture import CircularAperture, CircularAnnulus
from astropy import units as u



def calculate_background_area(data, src_x, src_y, optimal_radius, show_fig=False, save_fig=False):
    """
    Defines the background area for a circular annulus of width 20. Excludes dead area between 
    detectors.
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

    if show_fig:
        fig, ax = plt.subplots()
        ax.imshow(data, norm=matplotlib.colors.LogNorm())
        plot_polygon(final_region, ax=ax, add_points=False)
        plt.show()
        if save_fig:
            plt.savefig("background_area.pdf", dpi=1000)
    return final_region.area


def calculate_source_area(data, src_x, src_y, optimal_radius, show_fig=False, save_fig=False):
    """
    Defines the source area for the main source in an observation. Excludes dead area between detectors.
    """
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

    if show_fig:
        fig, ax = plt.subplots()
        ax.imshow(data, norm=matplotlib.colors.LogNorm())
        plot_polygon(final_region, ax=ax, add_points=False)
        aper.plot(ax)
        plt.show()
        if save_fig:
            plt.savefig("source_area.pdf", dpi=1000)
    return final_region.area

def source_counts(data, src_x, src_y, optimal_radius, show_fig=False, save_fig=False):
    """
    Performs photometry on the corresponding source area to get a net source count value.
    """
    aper = CircularAperture([src_x, src_y], optimal_radius)
    src_counts, cts_err = aper.do_photometry(data)
    
    if show_fig:
        fig, ax = plt.subplots()
        ax.imshow(data, norm=matplotlib.colors.LogNorm())
        aper.plot(ax)
        plt.show()
        plt.close()
        if save_fig:
            plt.savefig("source_counts.pdf", dpi=1000)
    return src_counts[0] * u.ct

def bkg_counts(data, src_x, src_y, optimal_radius, show_fig=False, save_fig=False):
    """
    Performs photometry on the corresponding background area to get a net background count value.
    """
    aper = CircularAnnulus([src_x, src_y], optimal_radius, optimal_radius + 20)
    bkg_counts, cts_err = aper.do_photometry(data)
    
    if show_fig:
        fig, ax = plt.subplots()
        ax.imshow(data, norm=matplotlib.colors.LogNorm())
        aper.plot(ax)
        plt.show()
        plt.close()
        if save_fig:
            plt.savefig("background_counts.pdf", dpi=1000)
    return bkg_counts[0] * u.ct