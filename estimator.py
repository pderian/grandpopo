#!/usr/bin/python

###########
#
# This is a demo script for cutyphoon "server edition".
# It allows to send data to/retrieve results from Typhoon through network sockets.
# Thus, it becomes possible to implement custom pre/post processing on data/motion.
#
# Please note:
#  - you need to be familiar with Typhoon's basic commands (-> userGuide.pdf, Section B.1)
#  - this mode is still [Work In Progress] and undocumented.
#  - use with care, no guarantee on this one.
#
# First, run the demo with:
#   python estimator.py -h
# to display usage information, parameters and stuff.
#
# Check main() function below for more details.
#
# Written by P. DERIAN 2014-09-15.
# (C) 2014-... P. DERIAN & CSU, Chico Research Foundation
# Do not distribute without permission.
#
###########


import datetime # time and stuff
import json #json I/O
import socket #network socket
import os #filesystem stuff
import re #regular expression
import sys #system stuff
###
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import numpy
import scipy.io as io
import scipy.interpolate as interpolate
###
from config import *

##################################
#### client/server parameters ####
##################################
#-Note: MUST NOT BE CHANGED
MESSAGESIZE = 32 #fixed length of messages, in bytes
MAXCHUNKSIZE = 4096 #fixed maximum chunk size for large data, in bytes
SEND_DTYPE = numpy.dtype('float32') #send data as float32
RECEIVE_DTYPE = numpy.dtype('float64') #receive as float64


####################################
###### COMMUNICATION FUNCTIONS #####
####################################

def send_array(array2D, sock, verb=False):

    """
    Send a 2D array through given socket.

    Arguments:
    array2D a 2D numpy array to be sent
    sock a valid socket
    verb enable verbose mode.

    Written by P. Derian 2013-05-22
    """

    arrayHeight = int(array2D.shape[0])
    arrayWidth = int(array2D.shape[1])
    arrayElementsize = array2D.dtype.itemsize
    arrayMemsize = arrayHeight*arrayWidth*arrayElementsize
    #size description
    sizeString = '{} {} {}'.format(arrayHeight, arrayWidth, arrayElementsize)
    #send size decriptor in fixed size string
    message = get_fixed_length(sizeString, MESSAGESIZE)
    sock.sendall(message)
    if verb:
        print "Sending {}x{} ({}b)".format(arrayHeight, arrayWidth, arrayMemsize)
    #send the raw binary data
    sock.sendall(array2D.tostring())

def receive_array(sock, dtype, verb=1):

    """
    Receive array through given socket.

    Arguments:
    dtype a numpy.dtype object for the incoming array
    sock a valid socket
    verb enable verbose mode.

    Return:
    the array

    Written by P. Derian 2013-05-22
    """

    #receive incoming size
    msg = sock.recv(MESSAGESIZE)
    if msg=='':
        raise RuntimeError("socket connection broken")
    if len(msg)!=MESSAGESIZE:
        print 'wrong message length ({}), expecting {}.'.format(len(msg), MESSAGESIZE)
        stopEstimation()
    elif verb>1:
        print " /received {}bytes".format(len(msg))
    dims = re.split('(\d+)', msg)
    #incoming data bytesize
    arrayHeight = int(dims[1])
    arrayWidth = int(dims[3])
    arrayElementsize = int(dims[5])
    #check data type conformity
    if arrayElementsize!=dtype.itemsize:
        print 'incoming data size is {} bytes, expecting {}'.format(arrayElementsize, dtype.itemsize)
        ServerProcess.closeall()
    #compute total size
    arrayMemsize = arrayHeight*arrayWidth*arrayElementsize
    if verb:
        print "Receiving {}x{} ({}b)".format(arrayHeight, arrayWidth, arrayMemsize)
    #collect data
    data = ''
    numbytesTotal = 0
    while(numbytesTotal<arrayMemsize):
            #how many bytes remain
            numbytes = arrayMemsize - numbytesTotal
            #correct w.r.t. to chunksize
            if numbytes>MAXCHUNKSIZE:
                numbytes = MAXCHUNKSIZE
            #receive at most numbytes
            msg = sock.recv(numbytes)
            #check validity
            if msg=='':
                raise RuntimeError("socket connection broken")
            #and actual received length (may differ frome expected)
            numbytes = len(msg)
            numbytesTotal += numbytes
            data += msg
            if verb>1:
                print " /received {}bytes, total={}".format(numbytes, numbytesTotal)

    return numpy.fromstring(data, dtype=dtype).reshape((arrayHeight, arrayWidth))

def receive_text(sock, verb=1):

    """
    Receive text through given socket.

    Arguments:
    dtype a numpy.dtype object for the incoming array
    sock a valid socket
    verb enable verbose mode.

    Return:
    the text

    Written by P. Derian 2014-08
    """

    #receive size
    msg = sock.recv(MESSAGESIZE)
    if msg=='':
        raise RuntimeError("socket connection broken")
    if len(msg)!=MESSAGESIZE:
        print 'wrong message length ({}), expecting {}.'.format(len(msg), MESSAGESIZE)
        stopEstimation()
    elif verb>1:
        print " /received {}bytes".format(len(msg))
    dims = re.split('(\d+)', msg)
    #incoming data size
    nChar = int(dims[1])
    if verb:
        print "Receiving {} chars".format(nChar)
    #collect data
    data = ''
    numbytesTotal = 0
    while(numbytesTotal<nChar):
            #how many bytes remain
            numbytes = nChar - numbytesTotal
            #correct w.r.t. to chunksize
            if numbytes>MAXCHUNKSIZE:
                numbytes = MAXCHUNKSIZE
            #receive at most numbytes
            msg = sock.recv(numbytes)
            #check validity
            if msg=='':
                raise RuntimeError("socket connection broken")
            #and actual received length (may differ frome expected)
            numbytes = len(msg)
            numbytesTotal += numbytes
            data += msg
            if verb>1:
                print " /received {}bytes, total={}".format(numbytes, numbytesTotal)

    return data

def receive_message(sock, verb=1):

    """
    Receive a message of MESSAGESIZE through given socket.

    Arguments:
    sock a valid socket
    verb enable verbose mode.

    Return:
    the message

    Written by P. Derian 2014-02-05.
    """
    #receive size
    msg = sock.recv(MESSAGESIZE)
    if msg=='':
        raise RuntimeError("socket connection broken")
    if len(msg)!=MESSAGESIZE:
        print 'wrong message length ({}), expecting {}.'.format(len(msg), MESSAGESIZE)
        sys.exit(-1)
    elif verb>1:
        print " /received {}bytes".format(len(msg))

    return msg

def get_fixed_length(message, fixedSize):

    """
    Transform given message into fixed length string (padding with spaces)

    Arguments:
    message the input message
    fixed size the desired size

    Written by P. Derian 2013-05-22
    """

    #check if message compatible with given size
    if len(message)>fixedSize:
        print 'message ({}) too long for fixed size {}.'.format(msg, fixedSize)
        sys.exit(-1)

    #descriptor format (including padding and fixed length)
    descrFormat='{{:<{}}}'.format(fixedSize)
    return descrFormat.format(message)

#################################
###### PREPROCESS FUNCTIONS #####
#################################

################################
###### ESTIMATION FUNCTION #####
################################
def estimate(infofile, serverPort, outputdir=None, comment='', verbose=True,
             save_probe=True, save_fields=True, coarse_resolution=1.):

    #-load data
    ###########
    if verbose:
        print 'Processing data from {}'.format(infofile)
    with open(infofile) as f:
        jsondata = json.load(f)
    # get its path for images
    # Note: this assumes images are in the same path as the infofile.
    datapath, _ = os.path.split(infofile)
    # some parameters
    width = jsondata['width']
    height = jsondata['height']
    n_frames = jsondata['numberFrame']
    tstart = datetime.datetime.strptime(jsondata['startTime'], '%Y-%m-%d %H:%M:%S')
    # the grid
    # Note: generate using the function in config.py
    # the dimensions should match 'width', 'height' attributes of jsondata.
    dx = jsondata['gridResolution']
    xx, yy = domain_grid(jsondata['gridOrigin'], jsondata['gridDimensions'],
                         jsondata['gridRotation'], jsondata['gridResolution'])
    if xx.shape[0]!=height or xx.shape[1]!=width:
        print 'Generated grid does not match image size! Abort.'
        sys.exit(-1)

    #-initialize timeseries data structure
    ######################################
    if save_probe:
        if verbose:
            print "Results will be save as timeseries in .json archive."
        # the probe
        is_probed = ((xx - AVG_PROBE['x'])**2 + (yy - AVG_PROBE['y'])**2) <= AVG_PROBE['r']**2
        # the data structure
        result_timeseries = {
            'method': 'wavelet-based optical flow (WOF) - (cu)Typhoon algorithm',
            'author': 'Pierre DERIAN',
            'website': 'www.pierrederian.net',
            'createdBy': __file__,
            'description': 'Apparent displacements in [pixel] estimated by (cu)Typhoon between images pairs contained in "sourceData". The image reference is used: origin is top-left, x positive towards right and y positive towards bottom. Displacements can be converted to (approximations of) instantaneous velocities by multiplying "ux" and "uy" by "dx"/"dt".',
            'comment': comment,
            'sourceData': os.path.basename(infofile),
            'label': jsondata['label'],
            'probe': AVG_PROBE,
            'grid': {'origin': jsondata['gridOrigin'],
                     'dimensions': jsondata['gridDimensions'],
                     'rotation': jsondata['gridRotation'],
                     'resolution': jsondata['gridResolution'],
                     },
            'ux': [],
            'uy': [],
            't': [],
            'dt': [],
            'dx': dx,
            }

    #-initialize fields data structure
    ##################################
    if save_fields:
        if verbose:
            print "Results will be save as fields in .mat archive."
        # the coarse grid
        dx_coarse = coarse_resolution
        if dx_coarse<dx:
            print "Decimated grid scale ({} m/px) should be greater than input image's ({} m/px). Abort.".format(dx_coarse, dx)
            sys.exit(-1)
        xx_coarse, yy_coarse = domain_grid(jsondata['gridOrigin'], jsondata['gridDimensions'],
                                           jsondata['gridRotation'], dx_coarse)
        height_coarse, width_coarse = xx_coarse.shape
        # the "regular grid" just for interpolators, in the image reference
        x_fine = dx*numpy.arange(width)
        y_fine = dx*numpy.arange(height)
        x_coarse = dx_coarse*numpy.arange(width_coarse)
        y_coarse = dx_coarse*numpy.arange(height_coarse)
        # the data structure
        result_fields = {
            'x': xx_coarse,
            'x_descr': 'UTM eastern (x) coordinates of grid points',
            'y': yy_coarse,
            'y_descr': 'UTM northern (y) coordinates of grid points',
            'vx': numpy.zeros((n_frames-1, height_coarse, width_coarse), dtype='float32'),
            'vx_descr': 'horizontal (x) velocity in [m/s], positive towards right (image reference)',
            'vy': numpy.zeros((n_frames-1, height_coarse, width_coarse), dtype='float32'),
            'vy_descr': 'vertical (y) velocity in [m/s], positive towards bottom (image reference)',
            'dx': dx_coarse,
            'dx_descr': 'spatial resolution in [meter/pixel]',
            'dt': [],
            'dt_descr': 'temporal resolution (time-step) in [s]',
            't': [],
            't_descr': 'datetime of each frame, format "%Y-%m-%d %H:%M:%S.%f". This is the time of the first frame of the pair used for each field estimation.',
            'sourceData': os.path.basename(infofile),
            'sourceData_descr': 'the json file containing the frames used for estimation',
            'label': jsondata['label'],
            'label_descr': 'a label for this experiment',
            'method': 'wavelet-based optical flow (WOF) - (cu)Typhoon algorithm',
            'author': 'Pierre Derian - contact@pierrederian.net',
            'createdBy': __file__,
            'description': 'Apparent velocities in [m/s] estimated by (cu)Typhoon between images pairs listed in "SourceData". For the motions, the image reference is used: origin is top-left, x positive towards right and y positive towards bottom. Fields were decimated from estimation resolution ({} m/px) to a coarser resolution ("dx"={} m/px) by interpolation.'.format(dx, dx_coarse),
            'comment': comment,
            }

    #-for each frame pair
    #####################
    for k in xrange(n_frames-1):
        # print image pair
        if verbose:
            print '\nEstimation {: 3d}/{}'.format(k+1, n_frames-1)

        #-pre-processing goes here
        ###########################
        if verbose:
            print "+++ Preprocessing..."
        # load source im0, im1
        im0_src = pyplot.imread(os.path.join(datapath, jsondata['sourceImages'][k]))
        im1_src = pyplot.imread(os.path.join(datapath, jsondata['sourceImages'][k+1]))
        # load filtered im0, im1
        im0_filt = pyplot.imread(os.path.join(datapath, jsondata['filteredImages'][k]))
        im1_filt = pyplot.imread(os.path.join(datapath, jsondata['filteredImages'][k+1]))
        # get the difference
        im0_preproc = im0_src.astype(SEND_DTYPE) - im0_filt.astype(SEND_DTYPE)
        im1_preproc = im1_src.astype(SEND_DTYPE) - im1_filt.astype(SEND_DTYPE)
        # [DEBUG] display
        if 0:
            fig0, (ax0, ax1) = pyplot.subplots(1,2)
            ax0.imshow(im0_preproc, cmap='gray', vmin=-48., vmax=48.,
                       extent=[x[0], x[-1], y[-1], y[0]])
            ax1.imshow(im1_preproc, cmap='gray', vmin=-48., vmax=48.,
                       extent=[x[0], x[-1], y[-1], y[0]])
            ax0.plot(xx[is_probed], yy[is_probed], '+b')
            ax1.plot(xx[is_probed], yy[is_probed], '+b')
            ax0.set_xlim(x[0], x[-1])
            ax0.set_ylim(y[-1], y[0])
            ax1.set_xlim(x[0], x[-1])
            ax1.set_ylim(y[-1], y[0])
            pyplot.show()

        #-perform estimation
        ####################
        try:
            if verbose:
                print "+++ Estimation..."
            #-send images to server
            #######################
            # create socket
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            # connect to server
            sock.connect(('localhost', serverPort))
            # send arrays (make sure the array is 2D and with appropriate SEND_DTYPE type!)
            send_array(im0_preproc, sock, verb=verbose)
            send_array(im1_preproc, sock, verb=verbose)

            #-retrieve estimation results
            #############################
            # receive horizontal displacement component
            Ux = receive_array(sock, RECEIVE_DTYPE, verb=verbose)
            # receive vertical displacement component
            Uy = receive_array(sock, RECEIVE_DTYPE, verb=verbose)
            # and a message about the execution
            if verbose:
                # fval, gnorm
                print 'Functional value; Gradient norm:'
                print '\t'+receive_message(sock, verb=verbose)
            # close
            sock.close()

            #-post-processing goes here
            ###########################
            # ....
            # typically, conversion displacement => velocity
            # and/or filtering, visualization, export to whatever format...
            if verbose:
                print "+++ Post processing..."
            # the time of img0 and the time step
            t_img0 = tstart + datetime.timedelta(seconds=jsondata['frameTimestamps'][k])
            dt_img = jsondata['frameTimestamps'][k+1] - jsondata['frameTimestamps'][k]
            if save_probe:
                # append estimate to result_timeseries
                result_timeseries['ux'].append(numpy.mean(Ux[is_probed]))
                result_timeseries['uy'].append(numpy.mean(Uy[is_probed]))
                # and the corresponding time (of img0) and time step (between image pair)
                result_timeseries['t'].append(t_img0.strftime('%Y-%m-%d %H:%M:%S.%f'))
                result_timeseries['dt'].append(dt_img)
            # and to result fields
            if save_fields:
                # interpolate motion to coarse grid and store in result
                # Note: the coordinates order is the same as numpy's dimensions
                # Note: and we convert to velocity using the original dx (NOT dx_coarse)
                interpolator = interpolate.RectBivariateSpline(y_fine, x_fine, Ux)
                result_fields['vx'][k] = (dx/dt_img)*interpolator(y_coarse, x_coarse, grid=True)
                interpolator = interpolate.RectBivariateSpline(y_fine, x_fine, Uy)
                result_fields['vy'][k] = (dx/dt_img)*interpolator(y_coarse, x_coarse, grid=True)
                # and the corresponding time (of img0) and time step (between image pair)
                result_fields['t'].append(t_img0.strftime('%Y-%m-%d %H:%M:%S.%f'))
                result_fields['dt'].append(dt_img)

            # [DEBUG] display
            if 0:
                dimY, dimX = Ux.shape
                xMotion = numpy.arange(dimX)
                yMotion = numpy.arange(dimY)
                qstep = 10
                fig = pyplot.figure(figsize=(1281./90., 720./90.))
                ax = pyplot.subplot(111, xlabel='x [px]', ylabel='y [px]',
                                    title='frame #{}'.format(k))
                pi = ax.imshow(im0_preproc, aspect='equal', cmap='copper')
                pq = ax.quiver(xMotion[8:-8:qstep], yMotion[8:-8:qstep],
                               Ux[8:-8:qstep, 8:-8:qstep], Uy[8:-8:qstep, 8:-8:qstep],
                               pivot='middle', color='white',
                               units='xy', angles='xy',
                               )
                pyplot.show()
#                 ax.set_xlim(200, 600)
#                 ax.set_ylim(200, 600)
#                 outfile = '../tmp/{:03d}.png'.format(k)
#                 fig.savefig(outfile, dpi=90.)
#                 if verbose:
#                     print 'saved {}'.format(outfile)
#                 pyplot.close(fig)
        except Exception as e:
            print e
            print '[!!] Estimation error - jumping to next pair.'
            pass

    #-write archive
    ###############
    if outputdir is not None:
        basename, _ = os.path.splitext(os.path.basename(infofile))
        if save_probe:
            outputfile = os.path.join(outputdir, '{}_{}probe.json'.format(
                jsondata['label'], jsondata['prefix']))
            with open(outputfile, 'w') as f:
                json.dump(result_timeseries, f, indent=0)
            if verbose:
                print 'Saved motion timeseries archive {}'.format(outputfile)
        if save_fields:
            outputfile = os.path.join(outputdir, '{}_{}fields.mat'.format(
                jsondata['label'], jsondata['prefix']))
            io.savemat(outputfile, result_fields)
            if verbose:
                print 'Saved motion fields archive {}'.format(outputfile)

##################################
###### POSTPROCESS FUNCTIONS #####
##################################

########################
###### MAIN SCRIPT #####
########################
def main(argv):
    #-those packages are used by the main script only
    ###
    import argparse #parse arguments
    import glob #find files
    import os #path stuff
    import pprint
    ###
    import scipy.ndimage as ndimage #advanced image processing stuff
    import scipy.io as io #to load/save MATLAB .mat
    import matplotlib.image as image #to read images
    import matplotlib.pyplot as pyplot #to plot stuff


    #-Default Parameters
    ####################
    #-description
    description = """
----------------------------------------------------
 This is a demo script for cutyphoon_server.
 It enables custom I/O and pre/post processing.
 Note that cutyphoon_server must be already running;
 please provide the port at which it is listening.
----------------------------------------------------"""
    epilog = """
If cutyphoon_server server is not running,
please start it with typical parameters, e.g.:
    cutyphoon_server -h -nM 7 -nD 8 -nE 8 \\
        -r 1 -a 1e-2 -port SOME_PORT > SOME_LOG 2>&1 &
where
    SOME_PORT is the listening port - typically pick something
              greater than 20000 and lesser than 256^2.
    SOME_LOG is the text file for the text output of cutyphoon_server.
"""

    #-set-up parser
    ###############
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("-f", "--file", dest="infofiles", type=str, nargs='*', default=[], help="json info files to be processed")
    parser.add_argument("serverPort", type=int, help="cutyphoon_server listening port")
    parser.add_argument("-s", "--stop", dest="stopServer", help="[EXPERIMENTAL] stop server when done", action="store_true")
    parser.add_argument("-q", "--quiet", dest="verbose", help="enable quiet mode", action="store_false")
    parser.add_argument("-o", "--output", default=None, dest="outputDir", help=".mat output directory")
    parser.add_argument("-c", "--comment", default='', dest="comment", help="optional comment to be included to the archive")
    parser.add_argument("--probe", dest="save_probe", action="store_true",
                        help="save probe timeseries (default)")
    parser.add_argument("--no-probe", dest="save_probe", action="store_false",
                        help="do NOT save probe timeseries")
    parser.add_argument("--fields", dest="save_fields", action="store_true",
                        help="save fields as .mat")
    parser.add_argument("--no-fields", dest="save_fields", action="store_false",
                        help="do NOT save fields as .mat (default)")
    parser.add_argument("-r", "--resolution", default=0.5, type=float, dest="coarse_res",
                        help="coarse resolution [m/px] of the output fields, with --fields.")

    parser.set_defaults(verbose=True, stopServer=False, save_probe=True, save_fields=False)
    args = parser.parse_args(argv)

    #-welcome
    if args.verbose:
        print "\n\t*** (cu)Typhoon-server control script - GrandPopo edition ***"

    #-check output dir
    if (args.outputDir is not None) and (not os.path.isdir(args.outputDir)):
        print 'Creating output directory:', args.outputDir
        os.makedirs(args.outputDirs, 0755)

    #-perform estimation for each file
    ##################################
    if len(args.infofiles)<1:
        raise RuntimeError("no data files were provided!")
    for infofile in args.infofiles:
        estimate(infofile, args.serverPort, outputdir=args.outputDir,
        comment=args.comment, verbose=args.verbose, save_probe=args.save_probe,
        save_fields=args.save_fields, coarse_resolution=args.coarse_res)

    #-stop server?
    ##############
    #-remember: server is still running, has to be stopped manually (e.g. kill -2 <PID>)
    #-an EXPERIMENTAL feature stops the server by sending 2 empty arrays.
    if args.stopServer:
        if args.verbose:
            print "Attempting to stop server. If nothing happens, hit Ctrl-C and kill it manually."
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect(('localhost', args.serverPort))
        send_array(numpy.empty((0,0), dtype=SEND_DTYPE), sock, verb=args.verbose)
        send_array(numpy.empty((0,0), dtype=SEND_DTYPE), sock, verb=args.verbose)
        sock.close()

    #-aaand we're done
    ##################
    if args.verbose:
        print "\n\t*** End of demo. Happy estimation! ***"



if __name__=="__main__":
    main(sys.argv[1:])
