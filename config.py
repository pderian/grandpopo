"""Configuration for the Grand Popo surface current analysis experiment.

Written by P. DERIAN 2017-01-06.
"""
### Path
ROOT_RAWDATA_DIR = '/Volumes/LaCie_Mac/pderian/data_GPP/Videos'
ROOT_PREPROC_DIR = '/Volumes/LaCie_Mac/pderian/data_GPP/Frames'
# this is the directory where ffmpeg & co are installed, since libav sucks ass.
FFMPEG_DIR = '/opt/local/bin' #'ffmpeg/ffmpeg-git-20160511-64bit-static'

### Case parameters
PARAMS_COMP60 = {
    # Parameters for the 60x60 m area covering the ADV/Aquapro sensors
    'label': 'timeseries_60m', #a label for the frame directory
    'image_format': 'jpg', #format of output images
    'median_length': 5., #length of median filter in [m]
    'resolution': 0.1, #[m/px]
    'x_bounds': (370295., 370355.), # x (easting) domain bounds in [m]
    'y_bounds': (694053., 694113.), # y (northing) domain bounds in [m]
    }
PARAMS_COMP30 = {
    # Parameters for the 30x30 m area covering the ADV/Aquapro sensors
    'label': 'timeseries_30m', #a label for the frame directory
    'image_format': 'jpg', #format of output images
    'median_length': 5., #length of median filter in [m]
    'resolution': 0.1, #[m/px]
    'x_bounds': (370310., 370340.), # x (easting) domain bounds in [m]
    'y_bounds': (694075.5, 694105.5), # y (northing) domain bounds in [m]
    }
PARAMS_RIP120 = {
    # Parameters for the 120x60 m area covering the dye release
    'label': 'flashrip_120m',
    'image_format': 'jpg', #format of output images
    'median_length': 5., #length of median filter in [m]
    'resolution': 0.2, #[m/px]
    'x_bounds': (370240., 370360.), # x (easting) domain bounds in [m]
    'y_bounds': (694045., 694105.), # y (northing) domain bounds in [m]
    }

### Averaging probe for time-series
AVG_PROBE = {
    'x': 370325.0, # x (easting) coordinate in [m]
    'y': 694098.0, # y (northing) coordinate in [m]
    'r': 1.0 # probe radius in [m]
}

### Projection
# this H was estimated from the mapping in the mat file.
# see preprocessor.estimate_H_from_mapping()
DEFAULT_H = [
    [2.4795542480487102e+03, 6.7539967752787433e+01,
     -2.2440068592504461e+05],
    [1.3228681378707056e+03, 3.5955050912544941e+01,
     -1.1965298848769323e+05],
    [3.5719691693753153e-03, 9.7256605585734505e-05,
     -3.2313166347066136e-01]
    ]