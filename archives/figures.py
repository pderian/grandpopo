###
import datetime
import itertools
import os.path
import cPickle
import sys
###
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import netCDF4
import numpy
import scipy.io as scio
import scipy.stats as stats
###

def figtracer(tracerfile):
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
    w = 39./6. #39 picas to inches
    h = 4.2 #inches
    dpi = 150. #dpi
    cmap = cm.get_cmap('magma')
    matplotlib.rcParams.update({'font.size': 10}) #change all default font size
    ####
    
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
        ax = pyplot.subplot(s, title='({}) frame #{:03d}'.format(l, n+data['k0']))
        # get the files, load image
        _, ifile, _ = data['files'][n]
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
    ax.set_ylabel('y (px)') 
    ax = axes[1] #top right
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax = axes[2] #bottom left
    ax.set_xlabel('x (px)')
    ax.set_ylabel('y (px)')
    ax = axes[3] #bottom right
    ax.set_xlabel('x (px)')
    ax.yaxis.set_ticklabels([])
    pyplot.subplots_adjust(left=0.09, bottom=0.09, right=0.98, top=0.95, wspace=0.1, hspace=0.12)
    fig.savefig('tracers.png', dpi=dpi)    
    pyplot.show()
    
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
        
def figpreproc(beachraw_file, beachpreproc_file, droneraw_file, droneB_file, droneBmed_file,
               beachtitle='', dronetitle=''):
    """
    Illustrates the preprocessing of shore and UAV data.
    Input and median filtered images.
    """
    #### plot parameters
    w = 39./6. #39 picas to inches
    h = 6. #inches
    dpi = 150. #dpi
    cmap = cm.get_cmap('magma')
    matplotlib.rcParams.update({'font.size': 10}) #change all default font size
    ####
    fig = pyplot.figure(figsize=(w, h))
    
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
    ax1 = fig.add_subplot(221, title='(a){}'.format(beachtitle), xlabel='x (m)', ylabel='y (m)')
    p1 = ax1.imshow(beachraw, extent=beachextent)
    ax2 = fig.add_subplot(222, title='(b){}'.format(beachtitle), xlabel='x (m)', ylabel='y (m)')
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
        
    ### drone
    drone_xlim = (700., 1000.)
    drone_ylim = (400., 100.)
    # plot images
    ax3 = fig.add_subplot(223, title='(c){}'.format(dronetitle), xlabel='x (px)', ylabel='y (px)')
    p3 = ax3.imshow(droneraw)
    ax4 = fig.add_subplot(224, title='(d){}'.format(dronetitle), xlabel='x (px)', ylabel='y (px)')
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
    
    #### finish    
    pyplot.subplots_adjust(left=.07, bottom=.08, right=.99, top=.95, hspace=.3, wspace=.32)
    fig.savefig('processing.png', dpi=dpi) 
    pyplot.show()

def figpreproc2(beachraw_file, beachradon_file, beachmedian_file,
                droneraw_file, droneradon_file, droneB_file, droneBmed_file,
                beachtitle='', dronetitle=''):
    """
    Illustrates the preprocessing of shore and UAV data.
    Input, radon and median filtered images.
    """
    #### plot parameters
    w = 39./6. #39 picas to inches
    h = 9. #inches
    dpi = 150. #dpi
    cmap = cm.get_cmap('magma')
    matplotlib.rcParams.update({'font.size': 10}) #change all default font size
    ####
    fig = pyplot.figure(figsize=(w, h))
    
    ### load data
    # beach
    beachraw = pyplot.imread(beachraw_file)
    beachradon_data = netCDF4.Dataset(beachradon_file)
    beachradon = beachradon_data['img'][26]
    beachmedian_data = netCDF4.Dataset(beachmedian_file)
    beachmedian = beachmedian_data['img'][26]
    beachextent = [0., beachmedian_data['x'][-1]-beachmedian_data['x'][0],
                   beachmedian_data['y'][-1]-beachmedian_data['y'][0], 0.] #in image reference                                  
    beacht0 = netCDF4.num2date(beachmedian_data.variables['t'][0],
                               beachmedian_data.variables['t'].units)
    beachtime = beacht0 + datetime.timedelta(seconds=float(beachmedian_data.variables['dt'][25]))
    print 'beach:', beachtime
    # drone
    droneraw = pyplot.imread(droneraw_file)
    droneradon = pyplot.imread(droneradon_file)
    droneB = pyplot.imread(droneB_file)
    droneBmed = pyplot.imread(droneBmed_file)
    dronemedian = droneB - droneBmed
    

    ### beach
    # plot images
    ax1 = fig.add_subplot(321, title='(a){}'.format(beachtitle), xlabel='x (m)', ylabel='y (m)')
    p1 = ax1.imshow(beachraw, extent=beachextent)
    ax3 = fig.add_subplot(323, title='(c){}'.format(beachtitle), xlabel='x (m)', ylabel='y (m)')
    p3 = ax3.imshow(beachradon, extent=beachextent, cmap='gray')
    p3.set_clim(-.3, .3)
    ax5 = fig.add_subplot(325, title='(e){}'.format(beachtitle), xlabel='x (m)', ylabel='y (m)')
    p5 = ax5.imshow(beachmedian, extent=beachextent, cmap='gray')
    p5.set_clim(-.3, .3)
    # style
    for ax in [ax1, ax3, ax5]:
        #ax.set_xlim(*drone_xlim)
        #ax.set_ylim(*drone_ylim)
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
        
    ### drone
    drone_xlim = (700., 1000.)
    drone_ylim = (400., 100.)
    # plot images
    ax2 = fig.add_subplot(322, title='(b){}'.format(dronetitle), xlabel='x (px)', ylabel='y (px)')
    p2 = ax2.imshow(droneraw)
    ax4 = fig.add_subplot(324, title='(d){}'.format(dronetitle), xlabel='x (px)', ylabel='y (px)')
    p4 = ax4.imshow(droneradon, cmap='gray', extent=drone_xlim+drone_ylim)
    #p4.set_clim(-.3, .3)
    ax6 = fig.add_subplot(326, title='(f){}'.format(dronetitle), xlabel='x (px)', ylabel='y (px)')
    p6 = ax6.imshow(dronemedian, cmap='gray')
    p6.set_clim(-.3, .3)
    # style
    for ax in [ax2, ax4, ax6]:
        ax.set_xlim(*drone_xlim)
        ax.set_ylim(*drone_ylim)
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(200))
    
    #### finish    
    pyplot.subplots_adjust(left=.07, bottom=.06, right=.99, top=.96, hspace=.3, wspace=.32)
    fig.savefig('processing_radon.png', dpi=dpi) 
    pyplot.show()

def figvectors():
    """
    Generates the vector plots
    """
    ### load data
    rootframedir = '/Users/pderian/Documents/Data/GrandPopo/plage/20140312'
    datafiles = ['/Users/pderian/Documents/Data/GrandPopo/plage/20140312/16.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/plage/20140312/17.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/plage/20140312/18.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/plage/20140312/19.mat',
                 '/Users/pderian/Documents/Data/GrandPopo/plage/20140312/20.mat',
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
    w = 39./6. #39 picas to inches
    h = 5. #inches
    dpi = 150. #dpi
    qstep = 20 #step for quiver
    matplotlib.rcParams.update({'font.size': 10}) #change all default font size
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
        axes.append(ax)
        # plot image
        p = ax.imshow(img, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], cmap='gray', zorder=0)
        # average fields over 10 s (~20 frames)
        #tmp_ux = ux[n,:,:]
        #tmp_uy = uy[n,:,:]
        ux_mean = numpy.mean(ux[n-10:n+10,:,:], axis=0)
        uy_mean = numpy.mean(uy[n-10:n+10,:,:], axis=0)
        tmp_ux = numpy.ma.array(ux_mean, mask=not_mask)
        tmp_uy = numpy.ma.array(uy_mean, mask=not_mask)
        # colors of the vectors
        vcolor = (YY<40)&(uy_mean<-0.75) 
        # plot vectors
        q = ax.quiver(x[::qstep], y[::qstep],
                      tmp_ux[::qstep, ::qstep], tmp_uy[::qstep, ::qstep],
                      vcolor[::qstep, ::qstep], cmap='RdYlBu_r',
                      units='xy', angles='xy', scale_units='xy', scale=.33,
                      width=.5, zorder=3)
        ax.set_xlim(27.5, 112.5)
        ax.invert_xaxis()
        ax.invert_yaxis()
    # tune individual axes
    ax = axes[0] # top left
    ax.plot(35., 45., 'xr', markeredgewidth=2., zorder=2)
    ax.set_ylabel('y (m)')
    ax = axes[2] # bottom left
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax = axes[3] # bottom right
    ax.set_xlabel('x (m)')
    # adjust layout
    pyplot.subplots_adjust(left=0.08, right=0.98, top=0.94, hspace=0.25)
    
    fig.savefig('dyevectors.png', dpi=dpi)
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
    
def figscatter(datafile):

    data = scio.loadmat(datafile)
    u_ADV = data['U_ADV'].squeeze()
    v_ADV = data['V_ADV'].squeeze()
    u_OF = data['U_OF'].squeeze()
    v_OF = data['V_OF'].squeeze()
    
    # filter out?
    tau = None
    if tau:
        good = (numpy.abs(u_OF)<tau) & (numpy.abs(v_OF)<tau) & (numpy.abs(u_ADV)<tau) & (numpy.abs(v_ADV)<tau) 
        u_ADV = u_ADV[good]
        v_ADV = v_ADV[good]
        u_OF = u_OF[good]
        v_OF = v_OF[good]
    
    # RMS error
    u_rms = numpy.sqrt(numpy.mean((u_ADV-u_OF)**2))
    v_rms = numpy.sqrt(numpy.mean((v_ADV-v_OF)**2))
    # linear regression
    u_slope, u_offset, u_r, u_p, u_stderr = stats.linregress(u_ADV, u_OF)
    v_slope, v_offset, v_r, v_p, v_stderr = stats.linregress(v_ADV, v_OF)
    
    print 'RMS: {:.2f}, {:.2f}'.format(u_rms, v_rms)
    print 'slopes: {:.2f}, {:.2f}'.format(u_slope, v_slope)
    print 'offsets: {:.2f}, {:.2f}'.format(u_offset, v_offset)
    print 'r^2: {:.2f}, {:.2f}'.format(u_r**2, v_r**2)
    print 'p: {:.2f}, {:.2f}'.format(u_p, v_p)
   #print 'stderr: {:.2f}, {:.2f}'.format(u_stderr, v_stderr) 
    
    fig = pyplot.figure()
    ax1 = fig.add_subplot(121, aspect='equal', xlabel='u ADV (m/s)', ylabel='u OF (m/s)')
    ax2 = fig.add_subplot(122, aspect='equal', xlabel='v ADV (m/s)', ylabel='v OF (m/s)')
    ax1.plot(u_ADV, u_OF, 'ob', markeredgecolor='none', alpha=0.05)
    ax2.plot(v_ADV, v_OF, 'or', markeredgecolor='none', alpha=0.05)
    for ax, sl, off, co in zip([ax1, ax2], [u_slope, v_slope], [u_offset, v_offset], ['b', 'r']):
        xlim = ax1.get_xlim()
        ylim = ax1.get_ylim()
        y = [sl*x + off for x in xlim]
        ax.plot(xlim, ylim, '--k')
        ax.plot(xlim, y, '--'+co)
    pyplot.show()

def figseries(datafile):
    data = scio.loadmat(datafile)
    for k, v in data.iteritems():
        if 'descr' in k:
            print k, v

if __name__=="__main__":

    # raw / filtered preprocess figure
    if 0:
        beachraw_file = '/Users/pderian/Documents/Data/GrandPopo/plage/20140313/1500_026.png'
        beachmedian_file = '/Users/pderian/Documents/Data/GrandPopo/plage/20140313/median_15_00.nc'
        beachradon_file = '/Users/pderian/Documents/Data/GrandPopo/plage/20140313/radon_15_00.nc'
        drone_frame = 115
        droneraw_file = "/Users/pderian/Documents/Data/GrandPopo/drone/halfres/rectif_jpg/{:03d}.jpg".format(drone_frame)
        droneradon_file = "/Users/pderian/Documents/Data/GrandPopo/drone/drone_{:03d}_radon.png".format(drone_frame)  
        droneB_file = "/Users/pderian/Documents/Data/GrandPopo/drone/halfres/rectif_B_png/{:03d}.png".format(drone_frame)  
        droneBmed_file = "/Users/pderian/Documents/Data/GrandPopo/drone/halfres/rectif_B_median/{:03d}.png".format(drone_frame)
        figpreproc2(beachraw_file, beachradon_file, beachmedian_file,
                   droneraw_file, droneradon_file, droneB_file, droneBmed_file)
                   #dronetitle='',# frame #{:03d}'.format(drone_frame),
                   #beachtitle='')
    
    # shore estimations - vector figures
    if 0:
        figvectors()
    if 0:
        shorecamTimestep()
        
    # shore estimations - ADV scatterplot
    if 1:
        datafile = '/Users/pderian/Documents/Data/GrandPopo/Data_Scatter_PLot_ADV_OF_GPP2.mat'
        figscatter(datafile)
    if 0:
        datafile = '/Users/pderian/Documents/Data/GrandPopo/series_median_5_20140313-mix.mat'
        figseries(datafile)
    
    # drone estimations - tracer figure  
    tracerfile = '/Users/pderian/Documents/tmp/gpp/drone/halfres_median/advect_avg_movie/tracerpath'
    if 0:
        figtracer(tracerfile)
    if 0:
        tracerlifetime(tracerfile)