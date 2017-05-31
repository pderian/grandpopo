"""Configuration for the Grand Popo surface current analysis experiment.

Written by P. DERIAN 2017-01-06.
"""
###
import numpy
###

### Path
ROOT_RAWDATA_DIR = '/Volumes/LaCie_Mac/pderian/data_GPP/Videos'
ROOT_PREPROC_DIR = '/Volumes/LaCie_Mac/pderian/data_GPP/Frames'
#ROOT_ESTIM_DIR = '/Volumes/LaCie_Mac/pderian/data_GPP/Estim/'
ROOT_ESTIM_DIR = '/Users/pderian/Documents/Data/GrandPopo/data/plage/estim'
# this is the directory where ffmpeg & co are installed, since libav sucks ass.
FFMPEG_DIR = '/opt/local/bin'

### The orientation of the beach reference w.r.t. Easting, Northing reference.
BEACH_ORIENTATION = 10. # in [degree]

### The local gravity at grand Popo
# source: https://www.sensorsone.com/local-gravity-calculator/
G_GRAV = 9.78089 # in [m/s2] at GrandPopo location (6.02 degree North)

### The water depth
H_WATER = 1. # in [m]

### Case parameters
PARAMS_COMP60 = {
    # Parameters for the 60x60 m area covering the ADV/Aquapro sensors (TGRS paper)
    'label': 'timeseries_60m', #a label for the frame directory
    'image_format': 'jpg', #format of output images
    'median_length': 5., #length of median filter in [m]
    'resolution': 0.1, #[m/px]
    'origin': (370295., 694053.), # origin (x,y) easting, northing of the domain in [m]
    'dimensions': (60., 60.), # domain size in [m]
    'rotation': 0., # domain rotation in [degree] around origin
    }
PARAMS_COMP30 = {
    # Parameters for the 30x30 m area covering the ADV/Aquapro sensors (TGRS paper)
    'label': 'timeseries_30m', #a label for the frame directory
    'image_format': 'jpg', #format of output images
    'median_length': 5., #length of median filter in [m]
    'resolution': 0.1, #[m/px]
    'origin': (370310., 694075.5), # origin (x,y) easting, northing of the domain in [m]
    'dimensions': (30., 30.), # domain size in [m]
    'rotation': 0., # domain rotation in [degree] around origin
    }
PARAMS_RIP120 = {
    # Parameters for the 120x60 m area covering the dye release (TGRS paper)
    'label': 'flashrip_120m',
    'image_format': 'jpg', #format of output images
    'median_length': 5., #length of median filter in [m]
    'resolution': 0.2, #[m/px]
    'origin': (370240., 694045.), # origin (x,y) easting, northing of the domain in [m]
    'dimensions': (120., 60.), # domain size in [m]
    'rotation': 0., # domain rotation in [degree] around origin
    }
PARAMS_SWASH125 = {
    # Parameters for the 125x45 m area covering most of the swash zone (Coastal Dyn paper)
    'label': 'swash_125m',
    'image_format': 'jpg', #format of output images
    'median_length': 5., #length of median filter in [m]
    'resolution': 0.2, #[m/px]
    'origin': (370240., 694050.), # origin (x,y) easting, northing of the domain in [m]
    'dimensions': (125., 45.), # domain size in [m]
    'rotation': BEACH_ORIENTATION, # domain rotation in [degree] around origin, see BEACH_ORIENTATION
    }
PARAMS_FIELD60 = {
    # Parameters for the 60x60 m area covering the ADV/Aquapro sensors (TGRS paper)
    # here showing unfiltered (color) images for the vector field movie
    'label': 'fields_60m', #a label for the frame directory
    'image_format': 'jpg', #format of output images
    'median_length': 0, #length of median filter in [m]
    'resolution': 0.1, #[m/px]
    'origin': (370295., 694053.), # origin (x,y) easting, northing of the domain in [m]
    'dimensions': (60., 60.), # domain size in [m]
    'rotation': 0., # domain rotation in [degree] around origin
    }

### Averaging probe for time-series
# Located near the instruments, between ADV and ADCP
# this probe was used in v0, v1 and v2.
AVG_PROBE_LEGACY = {
    'x': 370325.0, # x (easting) coordinate in [m]
    'y': 694098.0, # y (northing) coordinate in [m]
    'r': 1.0 # probe radius in [m]
}
# Centered on ADV
# this is used for the (vector field + time-series) video
AVG_PROBE_ADV = {
    'x': 370323.5, # x (easting) coordinate in [m]
    'y': 694097.8, # y (northing) coordinate in [m]
    'r': 1. # probe radius in [m]
}
# Set the default probe
AVG_PROBE = AVG_PROBE_LEGACY

### Projection
# this H was estimated from the mapping matrices defined in the mat file.
# see preprocessor.estimate_H_from_mapping()
DEFAULT_H = [
    [2.4795542480487102e+03, 6.7539967752787433e+01,
     -2.2440068592504461e+05],
    [1.3228681378707056e+03, 3.5955050912544941e+01,
     -1.1965298848769323e+05],
    [3.5719691693753153e-03, 9.7256605585734505e-05,
     -3.2313166347066136e-01]
    ]

### Helpers
def domain_grid(origin, dimensions, rotation, resolution):
    """Reference function for the generation of the domain grid from the origin, dimensions
    rotation and resolution parameters.

    Written by P. DERIAN 2017-01-11.
    """
    x = numpy.arange(0., dimensions[0], step=resolution,) # 1d
    y = numpy.arange(0., dimensions[1], step=resolution,) # 1d
    X, Y = numpy.meshgrid(x, y) # 2d
    cos_theta = numpy.cos(numpy.deg2rad(rotation))
    sin_theta = numpy.sin(numpy.deg2rad(rotation))
    Xr = cos_theta*X - sin_theta*Y + origin[0] #apply rotation
    Yr = sin_theta*X + cos_theta*Y + origin[1]
    return Xr, Yr
