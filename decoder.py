
###
import argparse
import datetime
import glob
import json
import os
import multiprocessing
import subprocess as sp
import sys
import time
###
import matplotlib.dates as dates
import matplotlib.pyplot as pyplot
import netCDF4
import numpy
import skimage.color as skcolor
import skimage.io as skio
###
from config import *
import preprocessor
###

#### FRAMES INPUT ####
######################
def video_info(videofile):
    """
    Retrieve width, height, number of frames and duration of the video file.
    See https://www.ffmpeg.org/ffprobe.html

    Written by P. DERIAN 2016-01-23
    """
    # command is like
    # ffprobe -select_streams v:0 -loglevel quiet -show_entries stream=index,width,height,nb_frames,duration -print_format json  myvideo.mpeg
    command = [os.path.join(FFMPEG_DIR,'ffprobe'),
               '-select_streams', 'v:0',
               '-loglevel', 'error',
               '-show_entries', 'format_tags=creation_time:stream=width,height,nb_frames,duration:frame=best_effort_timestamp_time',
               '-print_format', 'json',
               videofile,
               ]
    # run command
    pipe = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    infos = json.loads(pipe.stdout.read())
    pipe.terminate()
    # select datetime patten
    # because somehow it does not show up the same on different platforms
    if len(infos['format']['tags']['creation_time'])==19:
        time_value = infos['format']['tags']['creation_time']
        time_pattern = '%Y-%m-%d %H:%M:%S'
    elif len(infos['format']['tags']['creation_time'])>19:
        time_value = infos['format']['tags']['creation_time']
        time_pattern = '%Y-%m-%dT%H:%M:%S.%fZ' #not sure whyt the 'T' and 'Z'
    else:
        print '"creation_time" value: {} does not match any known pattern.'.format(
            infos['format']['tags']['creation_time'])
        sys.exit(-1)
    # finally return info
    return {'file': videofile,
            'width': int(infos['streams'][0]['width']),
            'height': int(infos['streams'][0]['height']),
            'nb_frames': int(infos['streams'][0]['nb_frames']),
            'duration': float(infos['streams'][0]['duration']),
            'creation_time': datetime.datetime.strptime(time_value, time_pattern),
            'timestamp': [float(f['best_effort_timestamp_time']) for f in infos['frames']],
            }

def load_keyframes(videofile, verbose=False):
    """
    Load all keyframes in the given video file.

    Arguments:
        video file: path to the video.
        verbose
    Output:
        a dict {'frames', 'info'} where
            frames is a sequence of frames
            info is a dict containing vido info.

    Written by P. DERIAN 2016-01-23
    """
    # Retrieve information on video content
    info = video_info(videofile)
    if verbose:
        print '\t{} frames ({}x{} px), {:.3f} s'.format(
            info['nb_frames'], info['width'], info['height'], info['duration']
            )
    nbytes = info['width']*info['height']*3
    # Extract frames
    # note: '-vsync 0' drops duplicates
    command = [os.path.join(FFMPEG_DIR,'ffmpeg'),
               '-loglevel', 'error',
               '-i', videofile,
               '-f', 'rawvideo',
               '-pix_fmt', 'rgb24',
               '-vsync', '0',
               'pipe:1',
               ]
    pipe = sp.Popen(command, stdout = sp.PIPE, bufsize=10**8)
    frames = []
    for k in xrange(info['nb_frames']):
        raw_image = pipe.stdout.read(nbytes)
        # transform the byte read into a numpy array
        image =  numpy.fromstring(raw_image, dtype='uint8')
        frames.append(image.reshape((info['height'],info['width'],3)))
    pipe.terminate()
    return {'frames':frames, 'info':info}

### PROCESSING FUNCTIONS ###
############################
def read_img_preprocess_save_img(inputfile, outputfile, preprocessor, verbose=True):
    """
    Preprocess the given image and save as a gray image

    Written by P. DERIAN 2016-01-23
    """
    # load and process
    skio.imsave(outputfile, preprocessor.grid_image(skio.imread(inputfile, as_grey=True)))

def read_video_preprocess_save_img(inputfile, outputdir, preprocessor, prefix='',
                                   verbose=True, workers=None, RGB=False, cmap='gray',
                                   format='png', label=''):
    """
    Read frames from the video file, preprocess and save as regular image(s).

    Written by P. DERIAN 2016-08-17
    Updated by P. DERIAN 2017-01-04: added json info output.
    Updated by P. DERIAN 2017-01-12: added grid rotation support.
    """
    ### load image data
    if verbose:
        print 'File', inputfile
        print '\tloading data...'
    data = load_keyframes(inputfile, verbose=verbose)

    ### process frames
    if verbose:
        print '\tprocessing...'
        sys.stdout.flush
    # allocate output
    processed_frames = []
    # create pool
    os.system("taskset -p 0xff %d" % os.getpid()) #see http://stackoverflow.com/questions/15639779/why-does-multiprocessing-use-only-a-single-core-after-i-import-numpy
    pool = multiprocessing.Pool(workers, maxtasksperchild=10)
    # create iterator
    def iterator(N, disp_progress=0):
        n = 0
        while n<N:
            # yield the RGB or a grayscale version of the frame
            # [TODO] enable to pass custom color=>grayscale converter?
            # Note: data is normalized to [0,1] here.
            yield ((data['frames'][n].astype(float)/255., True) if RGB
                   else (skcolor.rgb2grey(data['frames'][n]), True)
                   )
            # output progress
            if disp_progress and (not (n+1)%disp_progress):
                print '{} {:6.2f}%  (iter {})'.format(time.ctime(), 100.*(n+1.)/N, n+1)
            n += 1
    # process
    tic = time.time()
    result = pool.imap(
        preprocessor,
        iterator(N=len(data['frames']), disp_progress=10),
        chunksize=1,
        )
    pool.close()
    pool.join()
    # retrieve result
    for n, r in enumerate(result):
        processed_frames.append(r)
    toc = time.time()
    print 'done - {:.2f} s'.format(toc-tic)

    ### write images
    srcfiles = []
    filtfiles = []
    for i, images in enumerate(processed_frames):
        img_rect, img_filt = images
        # save the rectified version
        basename = 'src_{}{:03d}.{}'.format(prefix, i, format)
        outputfile = os.path.join(outputdir, basename)
        skio.imsave(outputfile, img_rect)
        print 'saved', outputfile
        srcfiles.append(basename)
        # save the filtered version, if any
        if img_filt is not None:
            basename = 'filt_{}{:03d}.{}'.format(prefix, i, format)
            outputfile = os.path.join(outputdir, basename)
            skio.imsave(outputfile, img_filt)
            print 'saved', outputfile
            filtfiles.append(basename)

    ### json data
    jsondata = {# metadata
                'author': 'Pierre DERIAN',
                'website': 'www.pierrederian.net',
                'createdBy': __file__,
                'description': 'Metadata for the frames extracted from video file and pre-processed. Resolution is in [m/px]',
                # source data
                'sourceVideo': os.path.basename(data['info']['file']),
                'numberFrame': data['info']['nb_frames'],
                'label': label,
                # output data
                'imageFormat': format,
                'sourceImages': srcfiles,
                'filteredImages': filtfiles,
                'prefix': prefix,
                # time info
                'startTime': data['info']['creation_time'].strftime('%Y-%m-%d %H:%M:%S'),
                'frameTimestamps': data['info']['timestamp'],
                # domain info
                'width': preprocessor.shape[1],
                'height': preprocessor.shape[0],
                'gridResolution': preprocessor.param['resolution'],
                'gridOrigin': preprocessor.param['origin'],
                'gridDimensions': preprocessor.param['dimensions'],
                'gridRotation': preprocessor.param['rotation'],
                # filters
                'medianLengthPx': preprocessor.param['median_length_px'],
                # projection matrix
                'projectionMatrix': preprocessor.H.tolist(),
                }
    jsonfile = os.path.join(outputdir, '{}info.json'.format(prefix))
    with open(jsonfile, 'w') as f:
        json.dump(jsondata, f, indent=0)
    print 'wrote {}'.format(jsonfile)

### MAIN ####
#############
def main(argv):
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default=None, dest="outputfile",
                        help="netcdf output file / image output directory (with -i, --img)")
    parser.add_argument("-m", "--median", default=0., type=float, dest="median_length",
                        help="length of median filter in [m]")
    parser.add_argument("-p", "--prefix", default='', dest="prefix",
                        help="output image prefix (with -i, --img)")
    parser.add_argument("-l", "--label", default='', dest="label",
                        help="output dataset label (with -i, --img)")
    parser.add_argument("-f", "--format", default='jpg', dest="format",
                        help="output image format (with -i, --img)")
    parser.add_argument("-O", "--origin", default=(370310., 694075.5), type=float, nargs=2,
                         dest="grid_origin", help="domain origin (x, y) in [m]")
    parser.add_argument("-d", "--dimensions", default=(30., 30.), type=float, nargs=2,
                         dest="grid_dimensions", help="domain dimensions (dim_x, dim_y) in [m]")
    parser.add_argument("-a", "--angle", default=0., type=float, dest="grid_rotation",
                        help="grid rotation around the origin in [degree]")
    parser.add_argument("-r", "--res", default=0.1, type=float, dest="grid_resolution",
                        help="grid resolution in [m/px]")
    parser.add_argument("videofile", type=str, help="video file to be processed")
    args = parser.parse_args(argv)

    # Create a preprocessor for images
    preproc = preprocessor.DataPreprocessor(
        H=preprocessor.DEFAULT_H,
        resolution=args.grid_resolution,
        origin=args.grid_origin,
        dimensions=args.grid_dimensions,
        rotation=args.grid_rotation,
        median_length=args.median_length,
        )

    # Process file
    read_video_preprocess_save_img(
        args.videofile,
        args.outputfile,
        preproc,
        prefix=args.prefix,
        format=args.format,
        label=args.label,
    )

if __name__=="__main__":
    main(sys.argv[1:])



