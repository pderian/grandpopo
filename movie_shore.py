"""Render frames for the shore-based surface current estimates video.

Written by P. DERIAN 2016-2017
www.pierrederian.net
"""
###
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
    rootframedir = '/Users/pderian/Documents/Data/GrandPopo/data/plage/20140312'
    datafiles = ['16.mat',
                 '17.mat',
                 '18.mat',
                 '19.mat',
                 '20.mat',
                 ]
    data = []
    framenames = []
    nfields = 0
    for f in datafiles:
        tmp = scio.loadmat(os.path.join(rootframedir, f), squeeze_me=True)
        print tmp.keys()
        data.append(tmp)
        nfields += tmp['dt'].size
        # generate what would be the input image name
        basename, _ = os.path.splitext(os.path.basename(str(tmp['file'])))
        framenames += ['{}_{:03d}.png'.format(basename, i) for i in xrange(tmp['dt'].size)]
    # "constant" var
    dx = tmp['dx']
    x = tmp['x']
    x -= x[0]
    y = tmp['y']
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
        ux[i,:,:] *= dx_over_dt[i] # displacements are now in [m/s]
        uy[i,:,:] *= dx_over_dt[i]
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
            ax = fig.add_axes([0.1, 0.05, .8, .9],
                              title='camera time: {:%Y-%m-%d %H:%M:%S.%f}'.format(t[n]),
                              xlabel='x (m)', ylabel='y (m)')
            # extra labels
            pyplot.figtext(0.5, 0.95, '{}'.format(title_prefix),
                           ha='center', va='top', fontsize='large')
            pyplot.figtext(0.9, 0.05, u'Pierre D\u00E9rian & Rafael Almar, 2017',
                           ha='right', va='bottom', fontsize='medium')
            # background frame
            img = pyplot.imread(os.path.join(rootframedir, framenames[n]))
            i = pyplot.imshow(img, origin='lower', extent=[x[0], x[-1], y[0], y[-1]])
            # draw only if dt "reasonable"
            pFlow.draw(ax, noTail=False, fadeTail=True, tailWidth=1.25)
            # reverse axes
            ax.invert_xaxis()
            ax.invert_yaxis()
            # save
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
        # multiply the advection velocity field by dt to get displacements in [m]
        pFlow.move(x, y, tmp_ux*dt[n], tmp_uy*dt[n], mask=mask)


def main():
    tracer('/Users/pderian/Documents/Data/GrandPopo/data/plage/v2_tmp',
           title_prefix='Surface current estimates by "Typhoon" optical flow\nFlash rip monitoring - Grand Popo, Benin\n')


if __name__=="__main__":
    main()
