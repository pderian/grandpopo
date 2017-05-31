"""Render frames for the vector + time-series movie.

Written by P. DERIAN 2017
www.pierrederian.net
"""
###
import datetime
import os
import sys
###
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
import matplotlib.dates as dates
import numpy
import pandas
import scipy.io as scio
import skimage.io as skio
import skimage.transform as sktransform
###
from config import *
from figures import load_ADV_timeseries
###
matplotlib.rcParams['text.usetex'] = False #because this was set to True in figures.py


def vector():

    def rotate(ux, uy, angle=-BEACH_ORIENTATION):
        cosp = numpy.cos(numpy.deg2rad(angle))
        sinp = numpy.sin(numpy.deg2rad(angle))
        vl = cosp*ux - sinp*uy
        vc = sinp*ux + cosp*uy
        return vl, vc

    #### plot parameters
    w = 1921. #[px]
    h = 1080. #[px]
    dpi = 90. #dpi
    fontsize = 12. #pt
    axcolor = '.2' #color of axes edges, ticks, labels, etc.
    color_ADV = 'orange'
    color_OFv = 'k'
    color_OFs = 'teal'
    # set the default font
    font = {'family' : 'sans-serif',
            'sans-serif':['Helvetica'],
            'weight' : 'normal',
            'size'   : fontsize}
    matplotlib.rc('font', **font)
    #
    qimin = 45 #area of displayed vector
    qimax = -15
    qjmin = 30
    qjmax = -30
    qstep = 2 #sub-sampling-step
    qscale = 1 #arrow scale
    tmin = datetime.datetime(2014,3,13,15,0,25) #bounds for time-series axes
    tmax = datetime.datetime(2014,3,13,15,3,25)

    ### image grid
    [ximg, yimg] = domain_grid(PARAMS_FIELD60["origin"], PARAMS_FIELD60["dimensions"],
                               PARAMS_FIELD60["rotation"], PARAMS_FIELD60["resolution"])

    ### the ADV
    # location
    projection = sktransform.ProjectiveTransform(numpy.array(DEFAULT_H))
    iyxADV = numpy.array([[360, 660],]) #pixels
    yxADV = projection(iyxADV) #meters
    # and data
    t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130413_refEastNorth.txt', rotate=0.)
    data_ADV = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)

    ### estimates
    imov = 0 #movie frame index
    ux_OF = [] #for the time-series
    uy_OF = []
    t_OF = []
    flist = ['/Volumes/LaCie_Mac/pderian/data_GPP/Estim/fields_60m/timeseries_60m_20140313_1500_fields.mat',
             '/Volumes/LaCie_Mac/pderian/data_GPP/Estim/fields_60m/timeseries_60m_20140313_1501_fields.mat',
             '/Volumes/LaCie_Mac/pderian/data_GPP/Estim/fields_60m/timeseries_60m_20140313_1502_fields.mat',
    ]
    for f in flist:
        # load estimated data
        data_OF = scio.loadmat(f, squeeze_me=True)
        yyyymmdd, hhmm, _ = str(data_OF['sourceData']).split('_') #dirty hack for image names
        yyyymmdd_hhmm = '{}_{}'.format(yyyymmdd, hhmm)
        img_basename = os.path.join('/Volumes/LaCie_Mac/pderian/data_GPP/Frames/fields_60m/20140313/20140313_15',
                                    yyyymmdd_hhmm, 'src_'+yyyymmdd_hhmm+'_{:03d}.jpg')
        # and the probe grid index
        is_probed = ((data_OF['x'] - AVG_PROBE_ADV['x'])**2 + (data_OF['y'] - AVG_PROBE_ADV['y'])**2) <= AVG_PROBE_ADV['r']**2

        for it in xrange(data_OF['t'].size):

            # the time
            tframe = datetime.datetime.strptime(data_OF['t'][it], '%Y-%m-%d %H:%M:%S.%f')
            tframe_next = tframe + datetime.timedelta(seconds=data_OF['dt'][it])
            # the ADV vector averaged over the inter-frame time-step
            tmp_ADV = data_ADV[tframe:tframe_next].mean()
            # the figure
            fig = pyplot.figure(figsize=(w/dpi,h/dpi))

            ### Left axes: vector plot
            ax = fig.add_axes([.07, .09, .45, .8],
                              xlabel='easting x (m)', ylabel='northing y (m)',
                              title='{0:%Y-%m-%d}\n{0:%H:%M:%S.%f} (frame #{1:03d})'.format(tframe, it))
            # the image
            imfile = img_basename.format(it)
            im = pyplot.imread(imfile)
            # the image
            p = ax.imshow(im, origin='lower',
                          extent=[ximg[0,0], ximg[0,-1],
                                  yimg[0,0], yimg[-1,0]],
                          zorder=1,
                          )
            # the probe
            ax.add_artist(
                patches.Circle((AVG_PROBE_ADV['x'], AVG_PROBE_ADV['y']), radius=AVG_PROBE_ADV['r'],
                               fill=True, color=color_OFs, alpha=0.66, zorder=4),
                )
            # the vectors
            q = ax.quiver(data_OF['x'][qimin:qimax:qstep, qjmin:qjmax:qstep],
                          data_OF['y'][qimin:qimax:qstep, qjmin:qjmax:qstep],
                          data_OF['vx'][it, qimin:qimax:qstep, qjmin:qjmax:qstep],
                          data_OF['vy'][it, qimin:qimax:qstep, qjmin:qjmax:qstep],
                          units='xy', angles='xy', pivot='middle',
                          scale_units='xy', scale=qscale, width=0.1,
                          color=color_OFv,
                          zorder=3,
                          )
            # the ADV
            qa = ax.quiver(yxADV[0,1], yxADV[0,0], tmp_ADV['vl'], tmp_ADV['vc'],
                           units='xy', angles='xy', pivot='middle',
                           scale_units='xy', scale=qscale, width=0.2,
                           color=color_ADV,
                           zorder=5,
                           )
            # the frame
            x0 = 370300. + 7.5 #origin
            y0 = 694060. + 10
            ux = numpy.array([5., 0.]) #reference vectors
            uy = numpy.array([0., 5.])
            vl, vc = rotate(ux, uy, angle=BEACH_ORIENTATION)
            q = ax.quiver([x0,x0], [y0,y0], vl, vc,
                          units='xy', angles='xy', pivot='tail',
                          scale_units='xy', scale=1, width=0.2,
                          color='w',
                          zorder=5,)
            ax.text(x0-1., y0+1., 'u', color='w', ha='center', va='center')
            ax.text(x0+1., y0-.75, 'v', color='w', ha='center', va='center')
            #
            ax.set_xlim(ximg[0,0]+10., ximg[0,-1]-10.)
            ax.set_ylim(yimg[0,0]+15., yimg[-1,0]-5.)
            # tune
            ax.invert_xaxis()
            ax.invert_yaxis()

            ### right axes: series
            axb = fig.add_axes([.62, .25, .35, .25],
                              xlabel='time (min:sec) after {:%H:%M}:00'.format(tmin),
                              ylabel='v (m/s)',
                              title='Along-shore')
            axt = fig.add_axes([.62, .65, .35, .25],
                              xlabel='time (min:sec) after {:%H:%M}:00'.format(tmin),
                              ylabel='u (m/s)',
                              title='Cross-shore')
            # ADV subseries, rotated
            tmp_ADV = data_ADV[tmin:tframe_next]
            vl_ADV, vc_ADV = rotate(tmp_ADV['vl'], tmp_ADV['vc'])
            axb.plot(tmp_ADV.index, vl_ADV, color=color_ADV, label='ADV (8 Hz)')
            axt.plot(tmp_ADV.index, vc_ADV, color=color_ADV)
            # OF subseries, rotated
            ux_OF.append(numpy.mean(data_OF['vx'][it][is_probed]))
            uy_OF.append(numpy.mean(data_OF['vy'][it][is_probed]))
            t_OF.append(tframe)
            tmp_OF = pandas.DataFrame({'vl':ux_OF, 'vc':uy_OF}, index=t_OF)
            vl_OF, vc_OF = rotate(tmp_OF['vl'], tmp_OF['vc'])
            axb.plot(tmp_OF.index, vl_OF, color=color_OFs, label='Typhoon OF (1-2 Hz)')
            axt.plot(tmp_OF.index, vc_OF, color=color_OFs)
            # tune axes
            for ax in [axb, axt]:
                ax.axhline(0., ls='--', color=axcolor, lw=0.75, zorder=10)
                ax.xaxis.set_major_locator(dates.SecondLocator([0, 30,]))
                ax.xaxis.set_minor_locator(dates.SecondLocator(range(0,60,10)))
                ax.xaxis.set_major_formatter(dates.DateFormatter('%M:%S'))
                ax.set_xlim(dates.date2num([tmin, tmax]))
                ax.set_ylim(-1.5, 1.5)

            ### labels
            axb.legend(loc='upper center', bbox_to_anchor=(0.5, -0.33),
                       frameon=False, ncol=2)
            pyplot.figtext(0.9, 0.04, u'Pierre D\u00E9rian & Rafael Almar, 2017',
                           ha='right', va='bottom', fontsize='medium')

            ### done
            fig.savefig('../TMP/{:03d}.png'.format(imov), dpi=dpi)
            imov += 1
            pyplot.close(fig)

if __name__=="__main__":
    vector()
