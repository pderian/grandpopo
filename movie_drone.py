"""Render the frames for the UAV video-based estimates.

Written by P. DERIAN 2016-2017
www.pierrederian.net
"""
###
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
import scipy.io as io
import scipy.interpolate as interpolate
import skimage.measure as measure
###
sys.path.append('/Users/pderian/Documents/Python/Tools')
import inr
import followTheFlow as flow
###
ROOTESTIMDIR = "/Users/pderian/Documents/Data/GrandPopo/data/drone/halfres/estim_median"
ROOTIMGDIR = "/Users/pderian/Documents/Data/GrandPopo/data/drone/halfres"

def edir(reldir):
    return os.path.join(ROOTESTIMDIR, reldir)
def idir(reldir):
    return os.path.join(ROOTIMGDIR, reldir)

def view(estimfile, imgfile=None, maskfile=None, outfile=None, step=10, title=''):
    # motion
    ux, uy = inr.readMotion(estimfile)
    x = numpy.arange(ux.shape[1])
    y = numpy.arange(uy.shape[0])
    # mask
    if maskfile is not None:
        mask = numpy.logical_not(pyplot.imread(maskfile).astype('bool'))
        ux = numpy.ma.array(ux, mask=mask)
        uy = numpy.ma.array(uy, mask=mask)

    dpi = 90.
    fig = pyplot.figure(figsize=(1920.5/dpi, 1080./dpi))
    if imgfile is not None:
        ax = fig.add_axes([0.1, 0.05, .8, .9], title=title)
        img = pyplot.imread(imgfile)
        #if img.ndim==3:
        #    img = numpy.average(img, axis=-1, weights=[.299, 0.587, 0.114])
        p = pyplot.imshow(img, cmap='magma')
        #p.set_clim(0., 1.)
    else:
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        norm = numpy.sqrt(ux**2 + uy**2)
        p = pyplot.imshow(norm, cmap='viridis')
        p.set_clim(0., 10.)
        pyplot.colorbar(p, shrink=.33)
        pyplot.subplots_adjust(left=0.05, right=0.98)

    q = ax.quiver(x[::step], y[::step], ux[::step, ::step], -uy[::step, ::step],
                      units='xy', color='k', pivot='middle', scale=1., width=1.)
    ax.set_xlim(200, 1100)
#     ax.set_ylim(1200, 600)

    if outfile is not None:
        fig.savefig(outfile, dpi=dpi)
        pyplot.close(fig)
    else:
        pyplot.show()

def tracer(estimlist, imglist=None, masklist=None, outdir=None, pColor=None, title_prefix='', k0=0):

    ### data
    Ne = len(estimlist)
    if imglist is None:
        imglist = Ne*[None,]
    if masklist is None:
        masklist = Ne*[None,]
    data = zip(estimlist, imglist, masklist)

    ### load first field
    efile, ifile, mfile = data[0]
    # motion
    ux, uy = inr.readMotion(efile)
    x = numpy.arange(ux.shape[1])
    y = numpy.arange(uy.shape[0])
    XX, YY = numpy.meshgrid(x, y)
    # image
    if ifile is not None:
        img = pyplot.imread(ifile)
    else:
        img = None
    # mask
    if mfile is not None:
        mask = pyplot.imread(mfile)
        # apply to field
        ux = ux*mask
        uy = uy*mask

    ### initial tracer field
    pMap = cm.ScalarMappable(norm=colors.Normalize(vmin=0., vmax=1., clip=True),
                             cmap='magma') if (pColor is None) else None
    #start box
    pNum = 100 #number of drifters
    seed = 19850131 #to ensure reproductible drifter paths
    pDomain = (200., 150., 1100., 370.) #advection domain
    pInit = (960., 200., 1000., 350.) #initialization domain ("seeding box")
    #initialize...
    pFlow = flow.ParticleFlow(pDomain, initBounds=pInit,
                              num=pNum, colormap=pMap, color=pColor,
                              maxLen=10, maxLife=float('inf'),
                              massExtinction=False, seed=seed, archive=True)

    ### for each time step
    k = k0
    for efile, ifile, mfile in data[1:]:
        if outdir is not None:
            # plot the current state
            dpi = 90.
            fig = pyplot.figure(figsize=(1920.5/dpi, 1080./dpi))
            ax = fig.add_axes([0.1, 0.05, .8, .9], title='{}frame #{:3d}'.format(title_prefix, k),
                              xlabel='m (px)', ylabel='n (px)')
            # extra labels
            pyplot.figtext(0.5, 0.98,
                           'Surface current estimates by "Typhoon" optical flow\nFlash rip monitoring by UAV - Grand Popo, Benin',
                           ha='center', va='top', fontsize='large')
            pyplot.figtext(0.9, 0.04, u'Pierre D\u00E9rian & Rafael Almar, 2017',
                           ha='right', va='bottom', fontsize='medium')
            if img is not None:
                i = pyplot.imshow(img, cmap='magma')
            ax.add_patch(patches.Rectangle((pDomain[0], pDomain[1]), pDomain[2]-pDomain[0], pDomain[3]-pDomain[1],
                                            edgecolor='red', facecolor='None', linewidth=0.5, linestyle='-'))
            ax.add_patch(patches.Rectangle((pInit[0], pInit[1]), pInit[2]-pInit[0], pInit[3]-pInit[1],
                                            edgecolor='blue', facecolor='None', linewidth=0.5))
            pFlow.draw(ax, noTail=True, fadeTail=True, tailWidth=1.5)
            ax.set_xlim(200, 1100)
            ax.set_ylim(500, 0)
            outfile = os.path.join(outdir, 'advect_{:03d}.png'.format(k-k0))
            fig.savefig(outfile, dpi=dpi)
            pyplot.close(fig)
            print 'saved', outfile
        # move particles
        pFlow.move(x, y, ux, uy, mask=mask, newInitBounds=pInit)
        # retrieve their average motion, shift the init box by the horizontal disp
        ux_box, uy_box = pFlow.last_median_displacement()
        pInit = (pInit[0]+ux_box, pInit[1]+0., pInit[2]+ux_box, pInit[3]+0.)
        # and to next frame
        k += 1
        ux, uy = inr.readMotion(efile)
        if ifile is not None:
            img = pyplot.imread(ifile)
        else:
            img = None
        # mask
        if mfile is not None:
            mask = pyplot.imread(mfile)
    ### now we have the archive
    archiveData = pFlow.archiveData
    archiveData['k0'] = k0
    archiveData['seed'] = seed
    archiveData['files'] = data
    archiveData['author']  = 'Pierre DERIAN - pierre.derian@inria.fr'
    archiveData['description'] = 'motion, image and mask files as well as tracer positions and colors.'
    archiveData['created_by'] = __file__
    outfile = os.path.join(outdir, 'tracerpath')
    with open(outfile, 'wb') as f:
        cPickle.dump(archiveData, f, protocol=2)
    print 'tracers archived in', outfile

def time_average(kmin=5, kmax=286, dk=5):
    """
        Rolling temporal average of motion (and masks).
    """
    def avg(flist, fout, mlist=None, mout=None):
        # motion
        ux = None
        uy = None
        for f in flist:
            ux_tmp, uy_tmp = inr.readMotion(f)
            ux = ux_tmp if (ux is None) else ux + ux_tmp
            uy = uy_tmp if (uy is None) else uy + uy_tmp
        ux /= len(flist)
        uy /= len(flist)
        inr.writeMotion(fout, ux, uy)
        # mask
        if (mlist is not None) and (mout is not None):
            mask = None
            for m in mlist:
                m_tmp = pyplot.imread(m) > 0.1
                mask = m_tmp if (mask is None) else mask & m_tmp
            mask = Image.fromarray(255*mask.astype('uint8'))
            mask.save(mout)
    # for each averaged frame
    for k in xrange(kmin, kmax):
        nmin = k - dk
        nmax = k + dk
        flist = [edir('{:03d}_UV.inr'.format(n)) for n in xrange(nmin, nmax+1)]
        mlist = [edir('{:03d}_Mask.png'.format(n)) for n in xrange(nmin, nmax+1)]
        fout = edir('average/{:03d}_avg_UV.inr'.format(k))
        mout = edir('average/{:03d}_avg_Mask.png'.format(k))
        avg(flist, fout, mlist, mout)


def main():
    action = 'advect'

    # display fields
    if action=='view':
        for k in xrange(5, 50):
            estim = edir('average/{:03d}_avg_UV.inr'.format(k))
            img = idir('rectif_jpg/{:03d}.jpg'.format(k))
            mask = edir('average/{:03d}_avg_Mask.png'.format(k))
            out = edir('quiv_avg/{}_{:03d}.png'.format(('quivimg' if img is not None else 'quivnorm'), k))
            view(estim, img, mask, outfile=out, title='frame #{:03d}'.format(k))
    # average fields
    elif action=='advect':
        imin = 5
        imax = 286
        # time-averaged fields
        estimlist = [edir('average/{:03d}_avg_UV.inr'.format(k)) for k in xrange(imin, imax)]
        masklist = [edir('average/{:03d}_avg_Mask.png'.format(k)) for k in xrange(imin, imax)]
        # instant fields
#         estimlist = [edir('{:03d}_UV.inr'.format(k)) for k in xrange(imin, imax)]
#         masklist = [edir('{:03d}_Mask.png'.format(k)) for k in xrange(imin, imax)]
        imglist = [idir('rectif_jpg/{:03d}.jpg'.format(k)) for k in xrange(imin, imax)]
        tracer(estimlist, imglist, masklist, outdir=idir('v2_tmp'), k0=imin,
               title_prefix='10-s averaged fields - ')


if __name__=="__main__":
    main()
