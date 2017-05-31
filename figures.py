"""Render figures for the TGRS paper.

Written by P. DERIAN 2017
www.pierrederian.net
"""
###
import datetime
import glob
import itertools
import json
import os.path
import cPickle
import sys
###
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.dates as dates
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import netCDF4
import numpy
import pandas
import skimage.io as skio
import scipy.io as scio
import scipy.signal as signal
import scipy.stats as stats
###
import preprocessor
from config import *

### HELPERS ###
###############
def set_axes_color(ax, color='k'):
    """Set axes elements (frame, ticks, labels, title, ...) to the given color.
    """
    ax.xaxis.label.set_color(color)
    ax.yaxis.label.set_color(color)
    ax.title.set_color(color)
    for spine in ax.spines.itervalues():
            spine.set_linewidth(0.75)
            spine.set_color(color)
    ax.tick_params(axis='both', which='both', color=color, labelcolor=color)

    return ax

def load_OF_timeseries(infofiles, rotate=0.):
    """Load the OF timesries from json files, apply reference rotation as requested.

    Arguments:
        - infofiles: the list of .json files;
        - rotate=0. an optioanl rotation angle in [degree].

    Returns: t, v_long, v_cross
        - t a sequence of timestamps;
        - v_long, v_cross: the long- and cross-shore velocities

    Written by P. DERIAN 2017-03.
    """
    if not len(infofiles):
        print 'load_OF_timeseries(): empty list of files'
        sys.exit(-1)
    # load data
    v = []
    t = []
    for infofile in infofiles:
        try:
            with open(infofile) as f:
                jsondata = json.load(f)
            # the datetimes
            tmp_t = [datetime.datetime.strptime(ts, '%Y-%m-%d %H:%M:%S.%f')
                     for ts in jsondata['t']]
            # the [px/frame] => [m/s] conversion factor
            tmp_factor = [jsondata['dx']/dt for dt in jsondata['dt']]
            # the velocity
            tmp_v = [[ux*f, uy*f] for ux,uy,f in zip(jsondata['ux'],
                                                     jsondata['uy'],
                                                     tmp_factor)]
            # now concatenate
            v += tmp_v
            t += tmp_t
        except Exception as e:
            print 'Loading {} failed with error:'.format(infofile), e
            print 'Moving on...'
            pass
    # as arrays
    v = numpy.array(v)
    # rotate the reference frame
    if rotate!=0.:
        cosp = numpy.cos(numpy.deg2rad(rotate))
        sinp = numpy.sin(numpy.deg2rad(rotate))
        v_long = cosp*v[:,0] - sinp*v[:,1]
        v_cross = sinp*v[:,0] + cosp*v[:,1]
    else:
        v_long = v[:,0]
        v_cross = v[:,1]
    return t, v_long, v_cross

def load_ADV_timeseries(datafile, rotate=0.):
    # load data
    datestr, timestr, ux, uy, _ = numpy.loadtxt(
        datafile,
        dtype=numpy.dtype('a10, a13, f8, f8, f8'),
        unpack=True,
        )
    # generate times...
    t = [datetime.datetime.strptime(tmp_d+' '+tmp_t, '%Y-%m-%d %H:%M:%S.%f')
         for tmp_d, tmp_t in zip(datestr, timestr)]
     # rotate the reference frame
    if rotate!=0.:
        cosp = numpy.cos(numpy.deg2rad(rotate))
        sinp = numpy.sin(numpy.deg2rad(rotate))
        v_long = cosp*ux - sinp*uy
        v_cross = sinp*ux + cosp*uy
    else:
        v_long = ux
        v_cross = uy
    return t, v_long, v_cross

def rolling_median_test(data, win='10S', tau=2.):
    # usual median test
    # discard if (residual w.r.t. median) > tau*(median of residuals)
    median = data.rolling(win).median()
    median['res'] = numpy.sqrt((data['vl']-median['vl'])**2 + (data['vc']-median['vc'])**2)
    median['mres'] = median['res'].rolling(win).median()
    isvalid = median['res']<tau*median['mres']
    return data[isvalid]

def rolling_median_test2(data, win='1T', tau=10):
    # variant median test
    # discard the tau % points with the highest residual w.r.t. median
    median = data.rolling(win).median()
    residual = numpy.sqrt((data['vl']-median['vl'])**2 + (data['vc']-median['vc'])**2)
    residual_cutoff = numpy.sort(residual)[-max(1, int(float(tau)*float(residual.size)/100.))]
    isvalid = residual < residual_cutoff
    return data[isvalid]

def welch_dopplershift(x, fs, nperseg):
    # 50% overlap
    noverlap = nperseg // 2
    # the segment indices
    step = nperseg - noverlap
    indices = numpy.arange(0, x.shape[-1]-nperseg+1, step)
    # for each segment
    phase_speed = numpy.sqrt(H_WATER*G_GRAV)
    result = []
    for ind in indices:
        # the segment data
        tmp_x = x[ind:ind+nperseg]
        # its periodgoram
        [tmp_f, tmp_Pxx] = signal.periodogram(tmp_x, fs, 'hanning', nperseg)
        # the Doppler shift
        tmp_mean = numpy.mean(tmp_x)
        tmp_df = tmp_f*(phase_speed/(tmp_x.mean() + phase_speed))
        result.append((tmp_f, tmp_df, tmp_Pxx))
    # regroup the various spectra
    fref = numpy.arange(0., fs, fs/nperseg)
    P0 = numpy.zeros_like(fref);
    nS = numpy.zeros_like(P0);
    for f, df, Pxx in result:
        ifds = numpy.digitize(df, fref)
        for j,k in enumerate(ifds):
            # if the destination bin (k) is valid
            if k>0 and k<P0.size:
                P0[k] += Pxx[j]
                nS[k] += 1.
    has_sample = nS>0
    fref = fref[has_sample]
    P0 = P0[has_sample]/nS[has_sample]
    return fref, P0, result

### MAIN PLOTTING FUNCTIONS ###
###############################
def figmap(as_grey=False, with_panels=False):

    def get_yx_polygon(X, Y):
        return numpy.array([[Y[iy, ix], X[iy, ix]]
                           for [iy, ix] in zip([0,-1,-1,0], [0,0,-1,-1])])

    #### plot parameters
    w = 7.16 #inches, IEEE double column
    h = 4.3 #inches
    dpi = 450. #dpi
    fontsize = 8. #pt
    axcolor = '.2' #color of axes edges, ticks, labels, etc.
    # set the default font
    font = {'family' : 'sans-serif',
            'sans-serif':['Helvetica'],
            'weight' : 'normal',
            'size'   : fontsize}
    matplotlib.rc('font', **font)

    #### preprocessor
    preprocess = preprocessor.DataPreprocessor(
        H=DEFAULT_H,
        origin=(370220., 694040.),
        dimensions=(150., 75.),
        rotation=0.,
        resolution=0.2,
        )

    ### elements
    # 60-m estimation area around sensors
    X60, Y60 = domain_grid(PARAMS_COMP60['origin'], PARAMS_COMP60['dimensions'],
                           PARAMS_COMP60['rotation'], PARAMS_COMP60['resolution'])
    yx60 = get_yx_polygon(X60, Y60)
    iyx60 = preprocess.projection.inverse(yx60) # pixel coords
    y60_label = yx60[:,0].min() # world coord for the text label
    x60_label = yx60[:,1].mean()
    area60_color = 'k'
    # 120x60 m wide area for flash rip
    X120, Y120 = domain_grid(PARAMS_RIP120['origin'], PARAMS_RIP120['dimensions'],
                             PARAMS_RIP120['rotation'], PARAMS_RIP120['resolution'])
    yx120 = get_yx_polygon(X120, Y120)
    iyx120 = preprocess.projection.inverse(yx120) # pixel coords
    y120_label = yx120[:,0].min() # world coord for the text label
    x120_label = yx120[:,1].mean()
    area120_color = 'k'
    # sensors and such
    iyxADV = numpy.array([[360, 660],]) #pixels
    yxADV = preprocess.projection(iyxADV)
    iyxRelease = numpy.array([[270, 1175],]) #pixels
    yxRelease = preprocess.projection(iyxRelease)
    # compas
    x0Compas = preprocess.X[0,0] + 110. #[m]
    y0Compas = preprocess.Y[0,0] + 80. #[m]
    lCompas = 10. #[m]
    x1Compas = x0Compas
    y1Compas = y0Compas + lCompas
    x2Compas = x0Compas + lCompas
    y2Compas = y0Compas
    yxCompas = numpy.array([[y0Compas, x0Compas],
                [y1Compas, x1Compas],
                [y2Compas, x2Compas]])
    iyxCompas = preprocess.projection.inverse(yxCompas)
    iyxCompas -= iyxCompas[0] - [[550, 1000]] #shift the arrows
    yxCompas -= yxCompas[0] - [[preprocess.Y[0,0]+5., preprocess.X[0,0]+5]]
    for i in [1,2]:
        iyxCompas[i,:] -= iyxCompas[0,:]
        iyxCompas[i,:] *= 30./numpy.sqrt(numpy.sum(iyxCompas[i,:]**2))

    #### load images
    framefile = 'resources/sample_frame_release.jpg'
    frame_img = skio.imread(framefile, as_grey=as_grey)

    ### figure
    fig = pyplot.figure(figsize=(w, h))
    # axes
    axfr = fig.add_axes([.08, .48, .66, .49], xlabel='$m$ (px)', ylabel='$n$ (px)')
    # style axes
    for ax in [axfr,]:
        set_axes_color(ax, axcolor)

    ### main frame
    # plot image
    axfr.imshow(frame_img, cmap='gray', interpolation='nearest', vmin=0.1, vmax=0.9)
    # plot domains
    axfr.add_artist(patches.Polygon(numpy.roll(iyx60,1,axis=-1), fill=False,
        color=area60_color, ls='-', lw=1.))
    axfr.add_artist(patches.Polygon(numpy.roll(iyx120,1,axis=-1), fill=False,
        color=area120_color, ls='--', lw=1.))
    # plot sensors
    axfr.plot(iyxADV[0,1], iyxADV[0,0], 'ow', markeredgewidth=.5)
    axfr.plot(iyxRelease[0,1], iyxRelease[0,0], 'dw', markeredgewidth=.5)
    # plot compas
    axfr.quiver([iyxCompas[0,1],],
               [iyxCompas[0,0],],
               [iyxCompas[1,1],],
               [iyxCompas[1,0],],
               color=['w', 'k'],
               units='xy', angles='xy',
               scale_units='xy',
               )
    axfr.text(iyxCompas[0,1]-5, iyxCompas[0,0]-5, 'N', color='w', weight='bold',
              va='bottom', ha='right')
    # misc
    axfr.set_xlim(0, frame_img.shape[1])
    axfr.set_ylim(frame_img.shape[0], 0)
    # reset formatter for ticks
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

    if with_panels:
        ### read images
        ADVfile = 'resources/grandpopo_ADV.jpg'
        dyefile = 'resources/grandpopo_dyerelease.jpg'
        copterfile = 'resources/grandpopo_hexacopter.jpg'
        mapfile = 'resources/grandpopo_map.jpg'
        camerafile = 'resources/grandpopo_camera2.jpg'
        ADV_img = skio.imread(ADVfile, as_grey=as_grey)
        dye_img = skio.imread(dyefile, as_grey=as_grey)
        copter_img = skio.imread(copterfile, as_grey=as_grey)
        map_img = skio.imread(mapfile, as_grey=as_grey)
        camera_img = skio.imread(camerafile, as_grey=as_grey)

        ### axes
        axADV = fig.add_axes([.02, .02, .2, .35], xlabel='', ylabel='')
        axdye = fig.add_axes([.27, .02, .33, .35], xlabel='', ylabel='')
        axcopter = fig.add_axes([.65, .02, .33, .35], xlabel='', ylabel='')
        axmap = fig.add_axes([.77, .675, .22, .3], xlabel='', ylabel='')
        axcamera = fig.add_axes([.77, .40, .22, .3], xlabel='', ylabel='')
        # style axes
        for ax in [axADV, axdye, axcopter, axmap, axcamera]:
            set_axes_color(ax, axcolor)

        ### main frame
        axfr.text(0.085, 0.95, "a)", color='w', ha='left', va='top',
                  transform=fig.transFigure, fontsize='large')

        ## ADV
        axADV.imshow(ADV_img, cmap='gray', interpolation='nearest', vmin=0.1, vmax=0.9)
        axADV.set_xticks([])
        axADV.set_yticks([])

        ## Dye
        axdye.imshow(dye_img, cmap='gray', interpolation='nearest', vmin=0.1, vmax=0.9)
        axdye.set_xticks([])
        axdye.set_yticks([])

        ## Hexacopter
        axcopter.imshow(copter_img, cmap='gray', interpolation='nearest', vmin=0.1, vmax=0.9)
        axcopter.set_xticks([])
        axcopter.set_yticks([])

        ## Map
        axmap.imshow(map_img, cmap='gray', interpolation='nearest', vmin=0.1, vmax=0.9)
        axmap.set_xticks([])
        axmap.set_yticks([])
        axmap.text(0.785, 0.95, "b)", color='w', ha='left', va='top',
                  transform=fig.transFigure, fontsize='large')

        ## camera
        axcamera.imshow(camera_img, cmap='gray', interpolation='nearest', vmin=0.1, vmax=0.9)
        axcamera.set_xticks([])
        axcamera.set_yticks([])
        axcamera.text(0.785, 0.62, "c)", color='k', ha='left', va='top',
                      transform=fig.transFigure, fontsize='large')

    fig.savefig('../figures/configuration_{:.0f}dpi.pdf'.format(dpi), dpi=dpi)
    fig.savefig('../figures/configuration_150dpi.png', dpi=150)

def figtracer(tracerfile, force_imgdir=None):
    """
    Plots the figure with virtual drifters superimposed on input UAV images
    """
    # plot all tracers
    def plot_tracers(ax, n=0, cmap=None):
        # for each tracer
        Nt = len(data['tracer_keys']) #number of tracers
        for i, p in enumerate(data['tracer_keys']):
            tracer = data[p]
            # the position
            pos = tracer['pos'][n]
            # the color in the movie, else color by id
            color = tracer['color'][n] if cmap is None else cmap(float(i)/(Nt-1.))
            ax.plot(pos[0], pos[1], '.', color=color, markersize=4)

    #### plot parameters
    w = 7.16 #inches, IEEE double column
    h = 4.2 #inches
    dpi = 450. #dpi
    fontsize = 8. #pt
    axcolor = '.2' #color of axes edges, ticks, labels, etc.
    cmap = cm.get_cmap('magma') #tracers color
    # set the default font
    font = {'family' : 'sans-serif',
            'sans-serif':['Helvetica'],
            'weight' : 'normal',
            'size'   : fontsize}
    matplotlib.rc('font', **font)

    # load data
    with open(tracerfile) as f:
        data = cPickle.load(f)
    # create figure
    fig = pyplot.figure(figsize=(w, h))
    axes = []
    # s is the subplot index
    # n the frame number in data
    # l the panel label
    for s, n, l in zip([221, 222, 223, 224],
                       [0, 45, 119, 235],
                       ['a', 'b', 'c', 'd']):
        # create axes
        ax = pyplot.subplot(s, title='({}) frame \#{:03d}'.format(l, n+data['k0']))
        set_axes_color(ax, axcolor) #change the elements color
        # get the files, load image
        _, ifile, _ = data['files'][n]
        if force_imgdir:
            ifile = os.path.join(force_imgdir, os.path.basename(ifile))
        im = pyplot.imread(ifile)
        # display image
        ax.imshow(im)
        # and tracers
        plot_tracers(ax, n, cmap=cmap)
        # style axes
        ax.set_xlim(580., 1020.) #pixels
        ax.set_ylim(380., 120.)
        # TODO:conversion pixel=>m???
        axes.append(ax)
    # specific styling
    ax = axes[0] #top left
    ax.xaxis.set_ticklabels([])
    ax.set_ylabel(r'$n$ (px)')
    ax = axes[1] #top right
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax = axes[2] #bottom left
    ax.set_xlabel(r'$m$ (px)')
    ax.set_ylabel(r'$n$ (px)')
    ax = axes[3] #bottom right
    ax.set_xlabel(r'$m$ (px)')
    ax.yaxis.set_ticklabels([])
    ### reset formatters
    for ax in axes:
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    pyplot.subplots_adjust(left=0.07, bottom=0.09, right=0.95, top=0.94, wspace=0., hspace=0.25)
    fig.savefig('../figures/tracers_drone_avg_{:.0f}dpi.png'.format(dpi), dpi=dpi)
    fig.savefig('../figures/tracers_drone_avg_150dpi.png', dpi=150)

def tracerlifetime(tracerfile):
    """
    Computes the lifetime of virtual tracers
    """
    # load data
    with open(tracerfile) as f:
        data = cPickle.load(f)
    # for each particle
    lifetimes = [];
    for k in data['tracer_keys']:
        particle = data[k]
        reset = [0,] + [i for i,c in enumerate(particle['color']) if c is None]
        reset.append(len(particle['color']))
        lifetimes += numpy.diff(reset).tolist()
    print numpy.mean(lifetimes), numpy.std(lifetimes)

def tracerpath(tracerfile, force_imgdir=None):

    ### load data
    with open(tracerfile) as f:
        data = cPickle.load(f)

    ### split the different paths
    paths = []
    # for each drifter
    for k in data['tracer_keys']:
        particle = data[k]
        # find when the particle was reset
        reset = [0,] + [i for i,c in enumerate(particle['color']) if c is None]
        reset.append(len(particle['color'])) # add the last time
        # and append the corresponding path to the list
        for k in xrange(len(reset)-1):
            tmp = numpy.array(particle['pos'][reset[k]:reset[k+1]])
            # only if it has at least 2 points
            if tmp.shape[0]>1:
                paths.append(tmp)

    ### plot
    # create axes
    ax = pyplot.gca()
    #set_axes_color(ax, axcolor) #change the elements color
    # get the files, load image
    _, ifile, _ = data['files'][119]
    if force_imgdir:
        ifile = os.path.join(force_imgdir, os.path.basename(ifile))
    im = pyplot.imread(ifile)
    # display image
    ax.imshow(im)
    # plot paths
    for p in paths:
        c = 'r' if p[-1,1]<p[0,1] else '.2'
        ax.plot(p[:,0], p[:,1], color=c)
    # style axes
    ax.set_xlim(580., 1020.) #pixels
    ax.set_ylim(380., 120.)
    pyplot.show()

def figpreproc(beachraw_file, beachpreproc_file, droneraw_file, droneB_file, droneBmed_file,
               beachtitle='', dronetitle=''):
    """
    Illustrates the preprocessing of shore and UAV data.
    Input and median filtered images.
    """
    #### plot parameters
    # double column
#     w = 7.16#in
#     h = 6. #inches
    # single column
    w = 3.5#in #39./6. #39 picas to inches
    h = 3.6#6. #inches
    dpi = 450. #dpi
    fontsize = 8. #pt
    axcolor = '.2' #color of axes edges, ticks, labels, etc.

    #### create figure
    # set the default font
    font = {'family' : 'sans-serif',
            'sans-serif':['Helvetica'],
            'weight' : 'normal',
            'size'   : fontsize}
    matplotlib.rc('font', **font)
    # create figure and axes
    fig = pyplot.figure(figsize=(w, h))
    ax1 = fig.add_subplot(221, title='(a){}'.format(beachtitle), xlabel=r'$x$ (m)', ylabel=r'$y$ (m)')
    ax2 = fig.add_subplot(222, title='(b){}'.format(beachtitle), xlabel=r'$x$ (m)', ylabel=r'$y$ (m)')
    ax3 = fig.add_subplot(223, title='(c){}'.format(dronetitle), xlabel=r'$x$ (px)', ylabel=r'$y$ (px)')
    ax4 = fig.add_subplot(224, title='(d){}'.format(dronetitle), xlabel=r'$x$ (px)', ylabel=r'$y$ (px)')
    for ax in [ax1, ax2, ax3, ax4]:
        set_axes_color(ax, axcolor)

    ### load data
    # beach
    beachraw = pyplot.imread(beachraw_file)
    beachpreproc_data = netCDF4.Dataset(beachpreproc_file)
    beachpreproc = beachpreproc_data['img'][26]
#     beachextent = [beachpreproc_data['x'][0], beachpreproc_data['x'][-1],
#                    beachpreproc_data['y'][-1], beachpreproc_data['y'][0]] #in UTM. COrrect orientation?
    beachextent = [0., beachpreproc_data['x'][-1]-beachpreproc_data['x'][0],
                   beachpreproc_data['y'][-1]-beachpreproc_data['y'][0], 0.] #in image reference
    beacht0 = netCDF4.num2date(beachpreproc_data.variables['t'][0],
                               beachpreproc_data.variables['t'].units)
    beachtime = beacht0 + datetime.timedelta(seconds=float(beachpreproc_data.variables['dt'][25]))
    print 'beach:', beachtime
    # drone
    droneraw = pyplot.imread(droneraw_file)
    droneB = pyplot.imread(droneB_file)
    droneBmed = pyplot.imread(droneBmed_file)
    dronepreproc = droneB - droneBmed

    ### beach
    # plot images
    p1 = ax1.imshow(beachraw, extent=beachextent)
    p2 = ax2.imshow(beachpreproc, extent=beachextent, cmap='gray')
    p2.set_clim(-.3, .3)
    # style
    for ax in [ax1, ax2]:
        #ax.set_xlim(*drone_xlim)
        #ax.set_ylim(*drone_ylim)
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
        ax.invert_xaxis() # to match with UTM orientation

    ### drone
    drone_xlim = (700., 1000.)
    drone_ylim = (400., 100.)
    # plot images
    p3 = ax3.imshow(droneraw)
    p4 = ax4.imshow(dronepreproc, cmap='gray')
    p4.set_clim(-.3, .3)
    # style
    for ax in [ax3, ax4]:
        ax.set_xlim(*drone_xlim)
        ax.set_ylim(*drone_ylim)
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(200))

    ### reset formatters
    for ax in [ax1, ax2, ax3, ax4]:
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

    #### finish
#     pyplot.subplots_adjust(left=.05, bottom=.08, right=.99, top=.95, hspace=.3, wspace=.24) #double column
    pyplot.subplots_adjust(left=.14, bottom=.08, right=.95, top=.96, hspace=.35, wspace=.5) #single column
    fig.savefig('../figures/processing_{:.0f}dpi.png'.format(dpi), dpi=dpi)
    fig.savefig('../figures/processing_150dpi.png', dpi=150)

def figvectors():
    """
    Generates the vector plots
    """
    ### load data
    rootframedir = '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140312'
    datafiles = ['/Users/pderian/Documents/Data/GrandPopo/data/plage/20140312/16.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140312/17.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140312/18.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140312/19.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140312/20.mat',
                 ]
    data = []
    framenames = []
    nfields = 0
    for f in datafiles:
        tmp = scio.loadmat(f)
        data.append(tmp)
        nfields += tmp['dt'].size
        # generate what would be the input image name
        basename, _ = os.path.splitext(os.path.basename(str(tmp['file'])))
        framenames += ['{}_{:03d}.png'.format(basename, i) for i in xrange(tmp['dt'].size)]
    # "constant" var
    dx = tmp['dx'].squeeze()
    x = tmp['x'].squeeze()
    x -= x[0]
    y = tmp['y'].squeeze()
    y -= y[0]
    XX, YY = numpy.meshgrid(x, y)
    # group fields and times (assuming same shape)
    ux = numpy.vstack(tuple([d['ux'] for d in data]))
    uy = numpy.vstack(tuple([d['uy'] for d in data]))
    dt = numpy.hstack(tuple([d['dt'].squeeze() for d in data]))
    t = list(itertools.chain.from_iterable([d['t'].tolist() for d in data]))
    # convert t to datetimes
    t = [datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S.%f') for s in t]
    # conversion factor
    dx_over_dt = dx/dt
    for i in xrange(len(dt)):
        ux[i,:,:] *= dt[i]
        uy[i,:,:] *= dt[i]
    # read mask
    mask = pyplot.imread(os.path.join(rootframedir, 'beach_sea_mask.png'))
    not_mask = numpy.logical_not(mask)

    ### [DEBUG] display frames # and times
#     for i, d in enumerate(t):
#         print i, '--', d
#     return

    ### display
    # plot parameters
    #w = 39./6. #39 picas to inches
    #h = 5. #inches
    w = 7.16 #inches
    h = 4.6 #inches
    dpi = 450. #dpi
    fontsize = 8. #pt
    axcolor = '.2' #color of axes edges, ticks, labels, etc.
    qstep = 20 #step for quiver
    # set the default font
    font = {'family' : 'sans-serif',
            'sans-serif':['Helvetica'],
            'weight' : 'normal',
            'size'   : fontsize}
    matplotlib.rc('font', **font)
    # custom colormap for vectors
    vec_cmap = colors.LinearSegmentedColormap.from_list(
        'vectors', [(0., 0., 0.), (1., 1., 1.)], N=2)

    # create figure
    fig = pyplot.figure(figsize=(w, h))
    axes = []
    # for each panel
    for s, n, l in zip([221, 222, 223, 224], #subplot number
                       [100, 180, 319, 420], #frame index - modified from Raf #120
                       #[115, 205, 355, 585],
                       ['(a)', '(b)', '(c)', '(d)'], #subplot title
                       ):
        print n, '/', t[n], framenames[n], (t[n+10]-t[n-10]).total_seconds()
        # load frame
        img = pyplot.imread(os.path.join(rootframedir, framenames[n]))
        # create subplot
        ax = fig.add_subplot(s, aspect='equal',
                             title='{} {:%H:%M:%S} UTC'.format(l, t[n]))
        set_axes_color(ax, axcolor) #manually change axes elements color
        axes.append(ax)
        # plot image
        p = ax.imshow(img, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], cmap='gray', zorder=0)
        # average fields over 10 s (~20 frames)
        #tmp_ux = ux[n,:,:] # these are the instantaneous fields
        #tmp_uy = uy[n,:,:]
        di = 10
        ux_mean = numpy.mean(ux[n-di:n+di,:,:], axis=0)
        uy_mean = numpy.mean(uy[n-di:n+di,:,:], axis=0)
        tmp_ux = numpy.ma.array(ux_mean, mask=not_mask)
        tmp_uy = numpy.ma.array(uy_mean, mask=not_mask)
        # colors of the vectors
        vcolor = (YY<40)&(uy_mean<-0.75)
        # plot vectors
        q = ax.quiver(x[::qstep], y[::qstep],
                      tmp_ux[::qstep, ::qstep], tmp_uy[::qstep, ::qstep],
                      vcolor[::qstep, ::qstep], cmap=vec_cmap,
                      units='xy', angles='xy', scale_units='xy', scale=.33,
                      width=.5, zorder=3)
        ax.set_xlim(27.5, 112.5)
        ax.set_ylim(10., 60.)
        ax.invert_xaxis()
        ax.invert_yaxis()
    # tune individual axes
    ax = axes[0] # top left
    ax.plot(35., 45., 'ow', markeredgewidth=.75, markeredgecolor=axcolor, zorder=2)
    ax.set_ylabel(r'$y$ (m)')
    ax = axes[2] # bottom left
    ax.set_xlabel(r'$x$ (m)')
    ax.set_ylabel(r'$y$ (m)')
    ax = axes[3] # bottom right
    ax.set_xlabel(r'$x$ (m)')
    # reset formatters
    for ax in axes:
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    # adjust layout
    pyplot.subplots_adjust(left=0.07, bottom=.08, right=0.98, top=0.94, hspace=0.25)
    fig.savefig('../figures/dyevectors_{:.0f}dpi.png'.format(dpi), dpi=dpi)
    fig.savefig('../figures/dyevectors_150dpi.png', dpi=150)
    pyplot.close(fig)

    # [WIP] display mean field
    fig = pyplot.figure(figsize=(w, h))
    # load frame
    img = pyplot.imread(os.path.join(rootframedir, framenames[319]))
    # create subplot
    ax = fig.add_subplot(111, aspect='equal',
                         title='Mean field')
    set_axes_color(ax, axcolor) #manually change axes elements color
    # plot image
    p = ax.imshow(img, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], cmap='gray', zorder=0)
    ux_mean = numpy.mean(ux, axis=0)
    uy_mean = numpy.mean(uy, axis=0)
    tmp_ux = numpy.ma.array(ux_mean, mask=not_mask)
    tmp_uy = numpy.ma.array(uy_mean, mask=not_mask)
    # plot vectors
    qstep = 10
    q = ax.quiver(x[::qstep], y[::qstep],
                  tmp_ux[::qstep, ::qstep], tmp_uy[::qstep, ::qstep],
                  color='w',
                  units='xy', angles='xy', scale_units='xy', scale=.33,
                  width=.5, zorder=3)
    ax.set_xlim(27.5, 112.5)
    ax.set_ylim(10., 60.)
    ax.invert_xaxis()
    ax.invert_yaxis()
    pyplot.show()

def shorecamTimestep():
    """
    Computes the distribution of the shore camera time-step.
    """
    ### load data
    rootframedir = '/Users/pderian/Documents/Data/GrandPopo/plage/20140312'
    datafiles = ['/Users/pderian/Documents/Data/GrandPopo/plage/20140312/16.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/plage/20140312/17.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/plage/20140312/18.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/plage/20140312/19.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/plage/20140312/20.mat',
                 ]
    dt = []
    nfields = 0
    for f in datafiles:
        tmp = scio.loadmat(f)
        dt.append(tmp['dt'].squeeze())
    dt = numpy.hstack(tuple(dt))
    mean_dt = numpy.mean(dt)
    std_dt = numpy.std(dt)
    fq = 1./dt
    mean_fq = numpy.mean(fq)
    std_fq = numpy.std(fq)
    print 'dt: mean={:.2f} s, std={:.2f} s - min={:.2f} s, max={:.2f} s'.format(
        mean_dt, std_dt, numpy.min(dt), numpy.max(dt))
    pyplot.hist(dt, bins=20, range=(0.25, 1.25), cumulative=True)
    pyplot.show()

def figseries(rotate=(-10.)):
    """

    Note: the difference (vc_surf - vc_below) with
        - vc_surf = vc_estim, the surface cross-shore velocity (by OF estim)
        - vc_bottom = vc_ADV, below-surface cross-shore velocity (by the ADV)
    remain +/- constant throughout the considered time window. Therefore the
    Doppler-corrected spectra are +/- the same whether using the overall average
    velocity, or by shifting each segment by the local average (welch_dopplershift).
    """

    #### plot parameters
    w1 = 3.5 #inches (43 picas) IEEE single colum
    w2 = 7.16 #inches (43 picas) IEEE double column
    h = 2. #inches
    dpi = 450. #dpi
    fontsize = 8. #pt
    axcolor = '.2' #color of axes edges, ticks, labels, etc.
    # set the default font
    font = {'family' : 'sans-serif',
            'sans-serif':['Helvetica'],
            'weight' : 'normal',
            'size'   : fontsize}
    matplotlib.rc('font', **font)

    ### load data
    print 'Vector rotation = {:.3f} (Reference rotation = {:.3f})'.format(rotate, -rotate)
    # load ADV
    print 'Loading ADV...'
    t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130413_refEastNorth.txt', rotate=rotate)
    tf_ADV = dates.date2num(t_ADV)
    # load estimates
    print 'Loading v2 estimates...'
    flist_30 = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_30m/timeseries_30m_20140313_1*_probe.json')))
    flist_60 = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_60m/timeseries_60m_20140313_1*_probe.json')))
    t_estim, vl_estim, vc_estim = load_OF_timeseries(flist_60, rotate=rotate)
    print '\t{} values loaded ({} - {})'.format(len(t_estim), t_estim[0], t_estim[-1])

    ### prepare series
    # make pandas objects
    data_ADV = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)
    data_estim = pandas.DataFrame({'vl':vl_estim, 'vc':vc_estim}, index=t_estim)
    # remove outliers
    #data_ADV = rolling_median_test(data_ADV, win='20S', tau=3.)
    #data_estim = rolling_median_test(data_estim, win='20S', tau=3.)

    # resample (and average), drop missing data
    resample_short = '2S' #2 s
    resample_large = '2T' #2 min
    data_ADV = data_ADV.resample(resample_short).mean().dropna()
    data_estim = data_estim.resample(resample_short).mean().dropna()
    # resample over larger intervals for correl
    data_ADV_l = data_ADV.resample(resample_large).mean().dropna()
    data_estim_l = data_estim.resample(resample_large).mean().dropna()
    # crop
    tmin = '2014-03-13 12:29'
    tmax = '2014-03-13 17:01'
    data_estim = data_estim[tmin:tmax]
    data_estim_l = data_estim_l[tmin:tmax]
    data_ADV = data_ADV[tmin:tmax]
    data_ADV_l = data_ADV_l[tmin:tmax]
    zmin = 1530 # this is for the zoom
    zmax = zmin + 300

    ### statistics
    tau = 0. #to filter out low velocities (0 =  no filter)
    correl_str = 'correl. coeff. r, r^2, p-value:'
    linreg_str = 'slope, offset:'
    err_str = 'RMSE, MRE:'
    for sample, d_ADV, d_estim in zip([resample_short, resample_large],
                                      [data_ADV, data_ADV_l],
                                      [data_estim, data_estim_l],
                                      ):
        print 'Data resampled at {}'.format(sample)
        for var, varlabel in zip(['vc', 'vl'],['u - cross-shore', 'v - along-shore']):
            # if we look at the 'instantaneous values'
            if (tau>0) and (sample==resample_short):
                # index of ADV values where var>tau
                isaboveTau = d_ADV[var][d_ADV[var].abs()>tau].index
                # intersection with estimates
                inboth = d_estim.index.intersection(isaboveTau)
                print '\t{} ({} points with |{}_ADV|>{}):'.format(varlabel, len(inboth), var, tau)
            else:
                inboth = d_estim.index.intersection(d_ADV.index)
                print '\t{} ({} points):'.format(varlabel, len(inboth), tau)
            # compute stats
            v_rmse = numpy.sqrt(numpy.mean((d_ADV[var][inboth] - d_estim[var][inboth])**2))
            v_mre = numpy.median(numpy.abs(d_ADV[var][inboth] - d_estim[var][inboth])/numpy.abs(d_ADV[var][inboth]))
            v_slope, v_offset, v_r, v_p, v_stderr = stats.linregress(d_ADV[var][inboth], d_estim[var][inboth])
            print '\t\t{:25s} {:.2f}, {:.2f}'.format(err_str, v_rmse, v_mre)
            print '\t\t{:25s} {:.2f}, {:.2f}, {:.4f}'.format(correl_str, v_r, v_r**2, v_p)
            print '\t\t{:25s} {:.2f}, {:.2f}'.format(linreg_str, v_slope, v_offset)

    ### plot spectrum
    # compute welch's estimates
    dt = 2. #must be the same as resample_short
    fs = 1./dt
    nwin = int(15*60*fs) #15 min x 60s/min x 1/dt sample/s
    offset = 10. #offset the cross-shore spectra
    _, P_vl_ADV = signal.welch(data_ADV['vl'], fs=fs, nperseg=nwin)
    _, P_vl_estim = signal.welch(data_estim['vl'], fs=fs, nperseg=nwin)
    f_ADV, P_vc_ADV = signal.welch(data_ADV['vc'], fs=fs, nperseg=nwin)
    f_estim, P_vc_estim = signal.welch(data_estim['vc'], fs=fs, nperseg=nwin)
    # compute global Doppler shift
    phase_speed = numpy.sqrt(H_WATER*G_GRAV) #phase speed in [m/s]
    mean_vc_estim = data_estim['vc'].mean()
    mean_vc_ADV = data_ADV['vc'].mean()
    mean_vl_estim = data_estim['vl'].mean()
    mean_vl_ADV = data_ADV['vl'].mean()
    f_vc_estim = f_estim*(phase_speed/(mean_vc_estim + phase_speed))
    f_vl_estim = f_estim*(phase_speed/(mean_vl_estim + phase_speed))
    f_vc_ADV = f_ADV*(phase_speed/(mean_vc_ADV + phase_speed))
    f_vl_ADV = f_ADV*(phase_speed/(mean_vl_ADV + phase_speed))
    print 'Mean vc: estim={:.3f}, ADV={:.3f} m/s'.format(mean_vc_estim, mean_vc_ADV)
    print 'Mean vl: estim={:.3f}, ADV={:.3f} m/s'.format(mean_vl_estim, mean_vl_ADV)
    # compute local Dopplershift
    f_ds_vc_estim, P_ds_vc_estim, _ = welch_dopplershift(data_estim['vc'], fs, nwin)
    f_ds_vl_estim, P_ds_vl_estim, _ = welch_dopplershift(data_estim['vl'], fs, nwin)
    f_ds_vc_ADV, P_ds_vc_ADV, _ = welch_dopplershift(data_ADV['vc'], fs, nwin)
    f_ds_vl_ADV, P_ds_vl_ADV, _ = welch_dopplershift(data_ADV['vl'], fs, nwin)
    # print stuff
    imaxOF = numpy.argmax(P_vc_estim)
    imaxADV = numpy.argmax(P_vc_ADV)
    print 'OF peak: f={:.3f}, period={:.3f}'.format(f_estim[imaxOF], 1./f_estim[imaxOF])
    print 'ADV peak: f={:.3f}, period={:.3f}'.format(f_ADV[imaxADV], 1./f_ADV[imaxADV])
    # create figure(s)
    #fig0 = pyplot.figure(figsize=(w1,h), dpi=90)
    #fig2 = pyplot.figure(figsize=(w1,h), dpi=90)
    fig3 = pyplot.figure(figsize=(w1,1.5*h), dpi=90)
    ax0 = fig3.add_subplot(211, xlabel='', ylabel=r'raw PSD (m\textsuperscript{2}/s)',
                          xscale='log', yscale='log')
    ax2 = fig3.add_subplot(212, xlabel=r'$f$ (Hz)', ylabel=r'Doppler-corrected PSD (m\textsuperscript{2}/s)',
                          xscale='log', yscale='log')
    set_axes_color(ax0, axcolor)
    set_axes_color(ax2, axcolor)
    # plot raw spectra
    ax0.plot(f_ADV[1:], P_vl_ADV[1:], '-', color='.4')
    ax0.plot(f_estim[1:], P_vl_estim[1:], '-', color='r')
    ax0.plot(f_ADV[1:], offset*P_vc_ADV[1:], '-', color='.4', lw=1.5)
    ax0.plot(f_estim[1:], offset*P_vc_estim[1:], '-', color='r', lw=1.5)
    # plot Doppler-shifted spectra
    if True:
        ax2.plot(f_vl_ADV[1:], P_vl_ADV[1:], '-', color='.4')
        ax2.plot(f_vl_estim[1:], P_vl_estim[1:], '-', color='r')
        ax2.plot(f_vc_ADV[1:], offset*P_vc_ADV[1:], '-', color='.4', lw=1.5)
        ax2.plot(f_vc_estim[1:], offset*P_vc_estim[1:], '-', color='r', lw=1.5)
    else:
        ax2.plot(f_ds_vl_ADV, P_ds_vl_ADV, '-', color='.4')
        ax2.plot(f_ds_vl_estim, P_ds_vl_estim, '-', color='r')
        ax2.plot(f_ds_vc_ADV, offset*P_ds_vc_ADV, '-', color='.4', lw=1.5)
        ax2.plot(f_ds_vc_estim, offset*P_ds_vc_estim, '-', color='r', lw=1.5)
    # tune limits
    ax0.set_xlim(1./900, 0.3)
    ax0.set_ylim(1e-1, 1.5e2)
    ax2.set_xlim(1./900, 0.3)
    ax2.set_ylim(1e-1, 1.5e2)
    # reset formatters
    ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    # tune margins
    #fig0.subplots_adjust(left=.17, bottom=.22, right=.99, top=.93)
    #fig2.subplots_adjust(left=.17, bottom=.22, right=.99, top=.93)
    fig3.subplots_adjust(left=.17, bottom=.12, right=.99, top=.98)
    # save original spectrum
    #outfile0 = '../figures/spectrum_rot{}{}_{:.0f}dpi.png'.format(
    #    ('p' if rotate>0 else 'm'), int(abs(rotate)), dpi)
    #fig0.savefig(outfile0,
    #             dpi=dpi)
    #print 'saved', outfile0
    #outfile0 = '../figures/spectrum_rot{}{}_150dpi.png'.format(
    #    ('p' if rotate>0 else 'm'), int(abs(rotate)))
    #fig0.savefig(outfile0,
    #             dpi=150)
    # save shifted spectrum
    #outfile2 = '../figures/spectrum_Dopplershift_rot{}{}_{:.0f}dpi.png'.format(
    #    ('p' if rotate>0 else 'm'), int(abs(rotate)), dpi)
    #fig2.savefig(outfile2,
    #             dpi=dpi)
    #print 'saved', outfile2
    #outfile2 = '../figures/spectrum_Dopplershift_rot{}{}_150dpi.png'.format(
    #    ('p' if rotate>0 else 'm'), int(abs(rotate)))
    #fig2.savefig(outfile2,
    #             dpi=150)
    # save both spectra
    outfile3 = '../figures/spectrum_duo_rot{}{}_{:.0f}dpi.png'.format(
        ('p' if rotate>0 else 'm'), int(abs(rotate)), dpi)
    fig3.savefig(outfile3,
                 dpi=dpi)
    print 'saved', outfile3
    outfile3 = '../figures/spectrum_duo_rot{}{}_150dpi.png'.format(
        ('p' if rotate>0 else 'm'), int(abs(rotate)))
    fig3.savefig(outfile3,
                 dpi=150)

    ### plot time-series
    for var, labeltop, labelbottom, ylim in zip(['vl', 'vc'],
                                                [r'along $\bar{v}$', r'cross $\bar{u}$'],
                                                [r'along $v$',  r'cross $u$'],
                                                [[-.3, 0.8], [-.3, .8]]):
        fig1 = pyplot.figure(figsize=(w2,h), dpi=90)
        ax1 = pyplot.subplot(211, ylabel='{} (m/s)'.format(labeltop),
                             title='',
                             )
        ax2 = pyplot.subplot(212, ylabel='{} (m/s)'.format(labelbottom),
                             xlabel='Video time - {:%Y-%m-%d} UTC'.format(data_estim.index[0]),
                             )
        # tune axes color
        for ax in [ax1, ax2]:
            set_axes_color(ax, axcolor)
        # plot stuff
        for data, color, label, lw, zorder in zip([data_estim, data_ADV],
                                                  ['r', '.4'], ['OF', 'ADV'],
                                                  [0.75, 0.75], [2, 1]):
            # rolling average over resample_large window (2 min)
            tmp1 = data.rolling(resample_large, center=False, min_periods=1).mean()
            # decimate for plot?
            #tmp1 = tmp1.resample('30S').first()
            # curves
            ax1.plot(tmp1.index, tmp1[var], color=color, label=label, lw=2*lw, zorder=5+zorder)
            ax2.plot(data.index, data[var], color=color, label=label, lw=lw, zorder=zorder)
        # labels
        for ax in [ax1, ax2]:
            ax.yaxis.set_label_coords(-0.05, 0.5)
        # x locators and formatters
        ax1.axhline(0., ls='--', color=axcolor, lw=0.75, zorder=10)
        ax2.axhline(0., ls='--', color=axcolor, lw=0.75, zorder=10)
        ax1.xaxis.set_major_locator(dates.MinuteLocator([0, 30]))
        ax1.xaxis.set_minor_locator(dates.MinuteLocator([10, 20, 40, 50]))
        ax1.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
        ax2.xaxis.set_major_locator(dates.MinuteLocator(range(0,60,2)))
        ax2.xaxis.set_minor_locator(dates.MinuteLocator(range(60)))
        ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
        ax1.set_xlim(dates.date2num(datetime.datetime.strptime(t, '%Y-%m-%d %H:%M'))
            for t in ['2014-03-13 12:31', '2014-03-13 16:59'])
        # y locators, formators
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
        ax1.yaxis.set_ticks([-0.5, 0., 0.5])
        ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
        ax2.yaxis.set_ticks([-1., 0., 1.])
        ### reset formatters
        for ax in [ax1, ax2]:
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        # limits
        ax1.set_ylim(ylim)
        ax1.axvspan(data_estim.index[zmin], data_estim.index[zmax], zorder=0,
                    facecolor='.8', edgecolor='none')
        ax2.set_xlim(data_estim.index[zmin], data_estim.index[zmax])
        ax2.set_ylim(-1.5, 1.5)
        # adjust
        pyplot.subplots_adjust(left=.07, bottom=.19, right=.97, top=.955, hspace=.26)
        outfile = '../figures/timeseries_{}_rot{}{}_{:.0f}dpi.png'.format(
            var, ('p' if rotate>0 else 'm'), int(abs(rotate)), dpi)
        fig1.savefig(outfile,
                     dpi=dpi)
        print 'saved', outfile
        outfile = '../figures/timeseries_{}_rot{}{}_150dpi.png'.format(
            var, ('p' if rotate>0 else 'm'), int(abs(rotate)))
        fig1.savefig(outfile,
                     dpi=150)

if __name__=="__main__":

    # domain / configuration
    if 0:
        figmap()

    # raw / filtered preprocess figure
    if 0:
        beachraw_file = '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140313/1500_026.png'
        beachmedian_file = '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140313/median_15_00.nc'
        beachradon_file = '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140313/radon_15_00.nc'
        drone_frame = 115
        droneraw_file = "/Users/pderian/Documents/Data/GrandPopo/data/drone/halfres/rectif_jpg/{:03d}.jpg".format(drone_frame)
        droneradon_file = "/Users/pderian/Documents/Data/GrandPopo/data/drone/drone_{:03d}_radon.png".format(drone_frame)
        droneB_file = "/Users/pderian/Documents/Data/GrandPopo/data/drone/halfres/rectif_B_png/{:03d}.png".format(drone_frame)
        droneBmed_file = "/Users/pderian/Documents/Data/GrandPopo/data/drone/halfres/rectif_B_median/{:03d}.png".format(drone_frame)
        figpreproc(beachraw_file, beachmedian_file,
                   droneraw_file, droneB_file, droneBmed_file)

    # shore estimations - vector figures
    if 0:
        figvectors()
    if 0:
        shorecamTimestep()

    # shore estimations - ADV scatterplot
    if 1:
        figseries()

    # drone estimations - tracer figure
    tracerfile = '/Users/pderian/Documents/Data/GrandPopo/data/drone/halfres/estim_median/advect_avg_movie/tracerpath'
    droneimg_dir = '/Users/pderian/Documents/Data/GrandPopo/data/drone/halfres/rectif_jpg'
    if 0:
        figtracer(tracerfile, force_imgdir=droneimg_dir)
    if 0:
        tracerlifetime(tracerfile)
    if 0:
        tracerpath(tracerfile, force_imgdir=droneimg_dir)
