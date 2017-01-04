
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
import preprocessor
###
# this is the time unit used in NetCDF files
DEFAULT_NC_TIME_UNIT = 'days since 2014-01-01 00:00:00'
# this is the directory where ffmpeg & co are installed, since libav sucks ass.
FFMPEG_DIR = '/opt/local/bin' #'ffmpeg/ffmpeg-git-20160511-64bit-static'

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
    pipe = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    infos = json.loads(pipe.stdout.read())
    pipe.terminate()
    return {'file': videofile,
            'width': int(infos['streams'][0]['width']),
            'height': int(infos['streams'][0]['height']),
            'nb_frames': int(infos['streams'][0]['nb_frames']),
            'duration': float(infos['streams'][0]['duration']),
            'creation_time': datetime.datetime.strptime(
                infos['format']['tags']['creation_time'], '%Y-%m-%d %H:%M:%S'),
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

### HELPERS ###
###############

def read_img_preprocess_save_img(inputfile, outputfile, preprocessor, verbose=True):
    """
    Preprocess the given image and save as a gray image

    Written by P. DERIAN 2016-01-23
    """
    # load and process
    skio.imsave(outputfile, preprocessor.grid_image(skio.imread(inputfile, as_grey=True)))

def read_video_preprocess_save_img(inputfile, outputdir, preprocessor, prefix='',
                                   verbose=True, workers=None, RGB=False, cmap='gray',
                                   format='jpg'):
    """
    Read frames from the video file, preprocess and save as a PNG.
    Note: this is a RGB version.

    Written by P. DERIAN 2016-08-17
    Updated by P. DERIAN 2017-01-04: added json info output.
    """
    ### load image data
    if verbose:
        print 'File', inputfile
        print '\tloading data...'
    data = load_keyframes(inputfile, verbose=verbose)

    ### json data
    jsondata = {# metadata
                'author': 'Pierre DERIAN',
                'website': 'www.pierrederian.net',
                'createdBy': __file__,
                'description': 'Metadata for the frames extracted from video file and pre-processed. Resolution is in [m/px]; xBounds and yBounds in [m] UTM.',
                # source data
                'sourceVideo': os.path.basename(data['info']['file']),
                'numberFrame': data['info']['nb_frames'],
                # time info
                'startTime': data['info']['creation_time'].strftime('%Y-%m-%d %H:%M:%S'),
                'frameTimestamps': data['info']['timestamp'],
                # domain info
                'resolution': preprocessor.param['resolution'],
                'xDimPx': preprocessor.x.size,
                'yDimPx': preprocessor.y.size,
                'xBoundsUTM': preprocessor.param['xbounds'],
                'yBoundsUTM': preprocessor.param['ybounds'],
                # filters
                'medianLengthPx': preprocessor.param['median_length_px'],
                # projection matrix
                'projectionMatrix': preprocessor.H.tolist(),
                }
    jsonfile = os.path.join(outputdir, '{}info.json'.format(prefix))
    with open(jsonfile, 'w') as f:
        json.dump(jsondata, f, indent=0)
    print 'wrote {}'.format(jsonfile)

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
            yield ((data['frames'][n],) if RGB else (skcolor.rgb2grey(data['frames'][n]),))
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
    for i, images in enumerate(processed_frames):
        img_rect, img_filt = images
        # save the rectified version
        outputfile = os.path.join(outputdir, 'src_{}{:03d}.{}'.format(prefix, i, format))
        pyplot.imsave(outputfile, img_rect, cmap=cmap)
        print 'saved', outputfile
        # save the filtered version, if any
        if img_filt is not None:
            outputfile = os.path.join(outputdir,
                                      'filt_{}{:03d}.{}'.format(prefix, i, format),
                                      )
            pyplot.imsave(outputfile, img_filt, cmap=cmap)
            print 'saved', outputfile

def read_video_preprocess_save_netcdf(inputfile, outputfile, preprocessor, verbose=True, workers=None):
    """
    """
    ### load image data
    if verbose:
        print 'File', inputfile
        print '\tloading data...'
    data = load_keyframes(inputfile, verbose=verbose)

    ### for each frame
    if verbose:
        print '\tprocessing...'
        sys.stdout.flush
    # allocate output
    src_frames = numpy.empty((len(data['frames']), preprocessor.y.size, preprocessor.x.size))
    filt_frames = numpy.empty((len(data['frames']), preprocessor.y.size, preprocessor.x.size))
    # create pool
    os.system("taskset -p 0xff %d" % os.getpid()) #see http://stackoverflow.com/questions/15639779/why-does-multiprocessing-use-only-a-single-core-after-i-import-numpy
    pool = multiprocessing.Pool(workers, maxtasksperchild=10)
    # create iterator
    def iterator(N, disp_progress=0):
        n = 0
        while n<N:
            # yield a greyscale version of the frame
            #yield (skcolor.rgb2grey(data['frames'][n]),)
            # yield the red layer
            yield (data['frames'][n][:,:,0]/255.,)
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
    # retrieve results
    for n, r in enumerate(result):
        src_img, filt_img = r
        processed_frames[n] = src_img - filt_img
    toc = time.time()
    print 'done - {:.2f} s'.format(toc-tic)

    ### create NetCDF file
    if outputfile is None:
        path, _ = os.path.splitext(inputfile)
        outputfile = path + '.nc'
    if verbose:
        print '\tcreating NetCDf', outputfile
    fdtype = numpy.dtype('float32') # this data type is used to convert to the desired precision

    ### root group
    # root contains all the main data
    rootgroup = netCDF4.Dataset(outputfile, 'w', format='NETCDF4')
    # create dimensions
    tDim = rootgroup.createDimension('t', 1)
    dtDim = rootgroup.createDimension('dt', processed_frames.shape[0])
    xDim = rootgroup.createDimension('x', processed_frames.shape[2])
    yDim = rootgroup.createDimension('y', processed_frames.shape[1])
    # create variables
    tVar = rootgroup.createVariable('t', 'f8', 't') #note this one has 'f8' (double precision)
    dtVar = rootgroup.createVariable('dt', 'f4', 'dt')
    xVar = rootgroup.createVariable('x', 'f4', 'x')
    yVar = rootgroup.createVariable('y', 'f4', 'y')
    imgVar = rootgroup.createVariable('img', 'f4', ('dt', 'y', 'x'))
    # set variable attributes
    tVar.units = DEFAULT_NC_TIME_UNIT
    tVar.description = 'date and time of the beginning of the video'
    dtVar.units = 'seconds'
    dtVar.description = 'timestamp of frames from the beginning of the video'
    xVar.units = 'meters'
    xVar.description = 'UTM eastern (x) coordinates of grid points'
    yVar.units = 'meters'
    yVar.description = 'UTM northern (y) coordinates of grid points'
    imgVar.units = 'luminance'
    imgVar.description = 'stack preprocessed and rectified data'
    # set global attributes
    rootgroup.author = 'Pierre Derian - contact@pierrederian.net'
    rootgroup.date_created = str(datetime.datetime.now())
    rootgroup.how_created = 'Created by {}'.format(__file__)
    rootgroup.description = 'This NetCDf archive contains data from the GrandPopo measurement campaign. Frames were extracted from the videofiles using "ffmpeg", preprocessed and rectified.',
    rootgroup.source_video = inputfile
    # fill data
    tVar[:] = netCDF4.date2num(data['info']['creation_time'], DEFAULT_NC_TIME_UNIT)
    dtVar[:] = numpy.array(data['info']['timestamp'], dtype=fdtype)
    xVar[:] = preprocessor.x.astype(fdtype)
    yVar[:] = preprocessor.y.astype(fdtype)
    imgVar[:] = processed_frames.astype(fdtype)
    ### preprocess group
    preprocgroup = rootgroup.createGroup('preprocess')
    # create dimensions
    hdim = preprocgroup.createDimension('h', 3)
    # create variables
    hVar = preprocgroup.createVariable('H', 'f8', ('h', 'h'))
    resVar = preprocgroup.createVariable('resolution', 'f4')
    radonVar = preprocgroup.createVariable('radon_length', 'f4')
    radonpxVar = preprocgroup.createVariable('radon_length_px', 'i4')
    medianVar = preprocgroup.createVariable('median_length', 'f4')
    medianpxVar = preprocgroup.createVariable('median_length_px', 'i4')
    # set attributes
    preprocgroup.description = 'Preprocessing parameters for the "img" data.'
    resVar.units = 'meters'
    resVar.description = 'resolution used for image rectification.'
    radonVar.units = 'meters'
    radonVar.description = 'physical length of the radon filter (0=disabled).'
    radonpxVar.units = 'pixels'
    radonpxVar.description = 'numerical length of the radon filter.'
    medianVar.units = 'meters'
    medianVar.description = 'physical length of the median filter (0=disabled).'
    medianpxVar.units = 'pixels'
    medianpxVar.description = 'numerical length of the median filter.'
    hVar.description = 'Matrix of the projective transform used for the rectification.'
    # set variables
    hVar[:] = preprocessor.H
    resVar[:] = preprocessor.param['resolution']
    radonVar[:] = preprocessor.param['radon_length']
    radonpxVar[:] = preprocessor.param['radon_length_px']
    medianVar[:] = preprocessor.param['median_length']
    medianpxVar[:] = preprocessor.param['median_length_px']
    ### end
    rootgroup.close()

### TEST ###
############
def test_image(preprocessor):
    """
    """
    files = ["/Volumes/Nymeria/Science/Data/GrandPopo/20140313/06/test/1.png",
             "/Volumes/Nymeria/Science/Data/GrandPopo/20140313/06/test/2.png",
             ]
    for i, f in enumerate(files):
        output = os.path.join('tmp', '{:03d}.png'.format(i))
        read_img_preprocess_save_img(f, output, preprocessor)

def test_video(preprocessor):
    """
    """
    files = ["/Volumes/Nymeria/Science/Data/GrandPopo/20140313/06/30.mp4",
             ]
    for f in files:
        data = load_keyframes(f, verbose=True)
        dt = data['info']['duration']/(data['info']['nb_frames']-1)
        for i, frame in enumerate(data['frames']):
            print i, i*dt
            #img = preprocessor.grid_image(skcolor.rgb2grey(frame))
            #skio.imsave(os.path.join('tmp', '{:03d}.png'.format(i)), img)
    return None

def test_frames_dt(preprocessor):

    file = sorted(glob.glob("/Volumes/Nymeria/Science/Data/GrandPopo/20140313/16/0*.mp4"))
    global_dt = []
    global_time = []
    global_info = []
    for f in file:
        info = video_info(f)
        global_info.append(info)
        t = [info['creation_time'] + datetime.timedelta(seconds=s) for s in info['timestamp']]
        global_time += t
        info['time'] = t
        print t[0], '->', t[-1]
        global_dt += [info['timestamp'][k+1]-info['timestamp'][k] for k in xrange(len(info['timestamp'])-1)]
    print global_dt[:15]
    nlarge = (numpy.asarray(global_dt)>0.6).sum()
    print 'dt > 0.6: {} / {} = {}%'.format(nlarge, len(global_dt), float(nlarge)/float(len(global_dt)))

    # histogram
    pyplot.figure()
    ax = pyplot.gca()
    bins = numpy.linspace(0, 1., num=21)
    ax.hist(global_dt, bins=bins)
    ax.set_xlabel('time between frame [s]')
    ax.set_ylabel('occurence')
    ax.set_title('GPP Video data - distribution of video time steps\n{:%Y-%m-%d %H:%M:%S} to {:%H:%M:%S}'.format(
        global_info[0]['time'][0], global_info[-1]['time'][-1]))
    # plot
    pyplot.figure()
    ax = pyplot.gca()
    ax.plot_date(global_time[:400], global_dt[:400], fmt='-b')
    ax.set_xlabel('time')
    ax.set_ylabel('time step [s]')
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    ax.set_title('GPP Video data - time steps\nstarting {:%Y-%m-%d %H:%M}'.format(
        global_info[0]['time'][0]))
    pyplot.show()

def test_netcdf(inputfile, outputdir, raw=False):
    """
    Check the content of a netcdf file of image data
    as created by read_video_preprocess_save_netcdf().

    Arguments:
        inputfile the input netcdf archive
        outputdir a directory for ouput images

    Written by P. DERIAN 2016-03-09
    """
    # load variables of interest
    rootgroup = netCDF4.Dataset(inputfile, 'r', format='NETCDF4')
    img = rootgroup.variables['img'][:]
    x = rootgroup.variables['x'][:]
    y = rootgroup.variables['y'][:]
    t0 = netCDF4.num2date(rootgroup.variables['t'][0], rootgroup.variables['t'].units)
    dt = rootgroup.variables['dt'][:]
    rootgroup.close()
    # display each frame
    for i in xrange(img.shape[0]):
        # restore frame datetime
        t = t0 + datetime.timedelta(seconds=float(dt[i]))
        if raw:
            outfile = outputdir+'_{:03d}.jpg'.format(i)
            pyplot.imsave(outfile, img[i], cmap='gray', vmin=0, vmax=1)
            print 'saved', outfile
        else:
            fig = pyplot.figure()
            ax = pyplot.axes(title=str(t), xlabel='UTM eastern [m]', ylabel='UTM northern [m]')
            ax.imshow(img[i], extent=[x[0], x[-1], y[-1], y[0]], interpolation='none', cmap='gray')
            fig.savefig(os.path.join(outputdir, '{:03d}.png'.format(i)))
            pyplot.close(fig)

### MAIN ####
#############
def main(argv):
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default=None, dest="outputfile", help="netcdf output file / image output directory (with -i, --img)")
    parser.add_argument("-r", "--radon", default=0., type=float, dest="radon_length", help="length of radon filter in [m]")
    parser.add_argument("-m", "--median", default=0., type=float, dest="median_length", help="length of median filter in [m]")
    parser.add_argument("--img", dest="asImg", action="store_true", help="export images")
    parser.add_argument("--nc", dest="asImg", action="store_false", help="export NetCDF (legacy)")
    parser.add_argument("-p", "--prefix", default='', dest="prefix", help="output image prefix (with -i, --img)")
    parser.add_argument("videofile", type=str, default=None, help="video file to be processed")
    parser.set_defaults(asImg=True)
    args = parser.parse_args(argv)

    # Parameters
    x_UTM = 370325 #sensors
    y_UTM = 694098
    #x_UTM = 370260 #flash rip
    #y_UTM = 694090
    dim = 30.
    resolution = .1 #sensors
    #resolution = .2 #flash rip
    radon_length = args.radon_length
    median_length = args.median_length
    xmin, xmax = x_UTM - dim/2., x_UTM + dim/2. #sensors
    #xmin, xmax = x_UTM - 20., x_UTM + 100. #flash rip
    ymin, ymax = y_UTM - 3*dim/4., y_UTM + dim/4.

    # Create a preprocessor for images
    preproc = preprocessor.DataPreprocessor(
        H=preprocessor.DEFAULT_H,
        xbounds=(xmin, xmax),
        ybounds=(ymin, ymax),
        resolution=resolution,
        radon_length=radon_length,
        median_length=median_length,
        )

    # Process file
    if args.asImg:
        read_video_preprocess_save_img(
            args.videofile,
            args.outputfile,
            preproc,
            prefix=args.prefix,
        )
    else:
        read_video_preprocess_save_netcdf(
            args.videofile,
            args.outputfile,
            preproc,
        )

    #p, _ =  os.path.splitext(args.outputfile)
    #test_netcdf(args.outputfile, p, True)

if __name__=="__main__":
    main(sys.argv[1:])



