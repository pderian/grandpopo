###
import argparse
import datetime
import glob
import json
import os
import sys
###
import matplotlib.pyplot as pyplot
import matplotlib.dates as dates
import numpy
import skimage.io as skio
###

def view(infofile, index=0):

    # load json data
    with open(infofile) as f:
        jsondata = json.load(f)
    # get its path for images
    path, _ = os.path.split(infofile)

    # load source image
    src_file = jsondata['sourceImages'][index]
    src_img = skio.imread(os.path.join(path, src_file))
    hp_img = src_img.astype(float)

    # load filtered image, if any
    if len(jsondata['filteredImages']):
        filt_file = jsondata['filteredImages'][index]
        filt_img = skio.imread(os.path.join(path, filt_file))
        hp_img -= filt_img.astype(float)
    else:
        filt_img = None
        hp_img = None

    # the grid
    x0, x1 = jsondata['xBoundsUTM']
    y0, y1 = jsondata['yBoundsUTM']
    dx = jsondata['resolution']
    x = dx*numpy.arange(jsondata['xDimPx'])
    y = dx*numpy.arange(jsondata['yDimPx'])
    #x = numpy.linspace(x0, x1, jsondata['xDimPx'])
    #y = numpy.linspace(y0, y1, jsondata['xDimPx'])

    # the time
    tstart = datetime.datetime.strptime(jsondata['startTime'], '%Y-%m-%d %H:%M:%S')
    timestamp = jsondata['frameTimestamps'][index]
    timg = tstart + datetime.timedelta(seconds=timestamp)

    # now display
    if filt_img is not None:
        ax1 = pyplot.subplot(131, title='rectified source image',
                             xlabel='x [m]', ylabel='y [m]')
        ax2 = pyplot.subplot(132, title='filtered image', xlabel='x [m]')
        ax3 = pyplot.subplot(133, title='high-pass image', xlabel='x [m]')
        p2 = ax2.imshow(filt_img, cmap='gray', vmin=0, vmax=255,
                        extent=[x[0], x[-1], y[-1], y[0]]) #note: reverse y order due to image reference
        p3 = ax3.imshow(hp_img, cmap='gray', vmin=-64, vmax=64,
                        extent=[x[0], x[-1], y[-1], y[0]])
    else:
        ax1 = pyplot.subplot(111, title='rectified source image')
    p1 = ax1.imshow(src_img, cmap='gray', vmin=0, vmax=255,
                    extent=[x[0], x[-1], y[-1], y[0]])
    pyplot.figtext(0.5, 0.98, '{} - frame #{:03d} - {}'.format(jsondata['sourceVideo'], index, timg),
                   va='top', ha='center')

    pyplot.show()

def load_timeseries(infofiles):

    # load data
    v = []
    t = []
    for infofile in infofiles:
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
    # as arrays
    v = numpy.array(v)
    t_plt = dates.date2num(t)
    # project to cross/longshore reference (-10 deg)
    cos10 = numpy.cos(numpy.deg2rad(10.))
    sin10 = numpy.sin(numpy.deg2rad(10.))
    v_long = cos10*v[:,0] + sin10*v[:,1]
    v_cross = -sin10*v[:,0] + cos10*v[:,1]

    # plot
#     figA = pyplot.figure()
#     ax1 = pyplot.subplot(211, ylim=(-4., 4.))
#     ax2 = pyplot.subplot(212, ylim=(-4., 4.))
#     ax1.axhline(0., ls='--', color='k', zorder=10)
#     ax2.axhline(0., ls='--', color='k', zorder=10)
#     ax1.plot_date(t_plt, v_long, '+', alpha=0.3)
#     ax2.plot_date(t_plt, v_cross, '+', alpha=0.3)
    #
#     figB = pyplot.figure()
#     pyplot.plot(v[:,0], v[:,1], '+', alpha=0.3)
    #
#     pyplot.show()
    return t, v_long, v_cross
    
def compare_timeseries():

    flist_30 = sorted(glob.glob('../data/plage/estim/timeseries_30m/timeseries_30m_20140313_15*_probe.json'))
    flist_60 = sorted(glob.glob('../data/plage/estim/timeseries_60m/timeseries_60m_20140313_15*_probe.json'))
    t_30, v_l30, v_c30 = load_timeseries(flist_30)
    t_60, v_l60, v_c60 = load_timeseries(flist_60)
    #
    fig = pyplot.figure()
    pyplot.subplot(121, aspect='equal')
    pyplot.plot(v_l30, v_l60, ',', alpha=0.1)
    pyplot.subplot(122, aspect='equal')
    pyplot.plot(v_c30, v_c60, ',', alpha=0.1)
    #
    hbins = numpy.linspace(-4, 4, num=161)
    fig = pyplot.figure()
    pyplot.subplot(121, aspect='equal')
    pyplot.hist2d(v_l30, v_l60, hbins, cmap='gray_r')
    pyplot.subplot(122, aspect='equal')
    pyplot.hist2d(v_c30, v_c60, hbins, cmap='gray_r')
    #
    pyplot.show()
    # [TODO] Pandas = > rolling mean. Shift time?
    # pour comparer à l'ADV, j'ai obtenu la meilleure correlation avec "tadv=tcam-377*0.5s"
    

if __name__=="__main__":
    #view('../tmp/test_v2/60m_jpg/20140613_0630_info.json')
    compare_timeseries()