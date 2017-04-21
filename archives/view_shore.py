import datetime
import itertools
import os.path
import cPickle
import sys
###
import PIL.Image as Image
import numpy
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
import scipy.io as scio
import scipy.interpolate as interpolate
import scipy.ndimage as ndimage
import skimage.measure as measure
###
sys.path.append('/Users/pderian/Documents/Python/Tools')
import inr
import followTheFlow as flow
###


def tracer(outdir=None, pColor=None, title_prefix='', average=False):

   
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
        ux[i,:,:] *= dt[i] # displacements are now in m/s
        uy[i,:,:] *= dt[i]
    # average in time
    ux_avg = ndimage.uniform_filter(ux, size=[20, 1, 1]) # about 10s centered window
    uy_avg = ndimage.uniform_filter(uy, size=[20, 1, 1])
    
    ### initial tracer field
    pMap = cm.ScalarMappable(norm=colors.Normalize(vmin=0., vmax=1., clip=False),
                             cmap='magma') if (pColor is None) else None
    pMap.cmap.set_over(alpha=0.) #hack to not display large motions
    #start box
    pNum = 1000
    seed = 19850131
    pDomain = (x[0], y[0], x[-1], y[-1])
    pInit = pDomain
    #initialize...
    pFlow = flow.ParticleFlow(pDomain, initBounds=pInit,
                              num=pNum, colormap=pMap, color=pColor,
                              maxLen=5, maxLife=20,
                              massExtinction=False, seed=seed, archive=False)
    # mask
    mask = pyplot.imread(os.path.join(rootframedir, 'beach_sea_mask.png'))

    ### for each time step
    for n in xrange(nfields):
        if outdir is not None:
            # plot the current state
            dpi = 90.
            fig = pyplot.figure(figsize=(1920.5/dpi, 1080./dpi))
            ax = fig.add_axes([0.1, 0.05, .8, .9], title='{}{:%Y-%m-%d %H:%M:%S.%f}'.format(title_prefix, t[n]),
                               xlabel='x (m)', ylabel = 'y (m)')
            img = pyplot.imread(os.path.join(rootframedir, framenames[n]))
            i = pyplot.imshow(img, origin='lower', extent=[x[0], x[-1], y[0], y[-1]])
            # draw only if dt "reasonable"
            pFlow.draw(ax, noTail=False, fadeTail=True)
#             ax.set_xlim(0, 120)
#             ax.set_ylim(0, 60)
            outfile = os.path.join(outdir, 'advect_{:03d}.png'.format(n))
            fig.savefig(outfile, dpi=dpi)
            pyplot.close(fig)
            print 'saved', outfile
        # move particles 
        if average:
            tmp_ux = ux_avg[n,:,:]
            tmp_uy = uy_avg[n,:,:]
        else:
            tmp_ux = ux[n,:,:]
            tmp_uy = uy[n,:,:]
        pFlow.move(x, y, tmp_ux*dt[n], tmp_uy*dt[n], mask=mask)


def main():
    tracer('/Users/pderian/Documents/Data/GrandPopo/plage/movie20140312_instant_r',
           title_prefix='"Typhoon" instantaneous estimates - Grand Popo, Benin - ')


if __name__=="__main__":
    main()