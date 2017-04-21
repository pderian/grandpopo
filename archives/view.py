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
import pandas
import scipy.interpolate as interpolate
import scipy.io as scio
import scipy.signal as signal
import scipy.stats as stats
import skimage.io as skio
###
from config import *

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

### time-series

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

def load_OF_timeseries(infofiles, rotate=0.):
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

def compare_timeseries(resample='1T'):

    flist_30 = sorted(glob.glob('../data/plage/estim/timeseries_30m/timeseries_30m_20140313_15*_probe.json'))
    flist_60 = sorted(glob.glob('../data/plage/estim/timeseries_60m/timeseries_60m_20140313_15*_probe.json'))
    flist_125 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/swash_125m/20140313/20140313_15/swash_125m_20140313_15*_probe.json'))

    # load data
    t_A, vl_A, vc_A = load_timeseries(flist_125, 10.)
    t_B, vl_B, vc_B = load_timeseries(flist_60, 0.)
    #dir_A = numpy.rad2deg(numpy.arctan2(vc_A, vl_A))
    #dir_B = numpy.rad2deg(numpy.arctan2(vc_B, vl_B))

    # RMS difference
    vl_rms = numpy.sqrt(numpy.mean((vl_A - vl_B)**2))
    vc_rms = numpy.sqrt(numpy.mean((vc_A - vc_B)**2))
    # linear regression
    vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(vl_A, vl_B)
    vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(vc_A, vc_B)
    print 'Raw instantaneous stats:'
    print '\tRMS: {:.2f}, {:.2f}'.format(vl_rms, vc_rms)
    print '\tslopes: {:.2f}, {:.2f}'.format(vl_slope, vc_slope)
    print '\toffsets: {:.2f}, {:.2f}'.format(vl_offset, vc_offset)
    print '\tr^2: {:.2f}, {:.2f}'.format(vl_r**2, vc_r**2)
    print '\tp: {:.2f}, {:.2f}'.format(vl_p, vc_p)

    #
    fig = pyplot.figure()
    ax1 = pyplot.subplot(121, aspect='equal')
    pyplot.plot(vl_A, vl_B, '.', alpha=0.1)
    xlim = ax1.get_xlim()
    pyplot.plot(xlim, xlim, '--k')
    ax2 = pyplot.subplot(122, aspect='equal')
    pyplot.plot(vc_A, vc_B, '.', alpha=0.1)
    xlim = ax2.get_xlim()
    pyplot.plot(xlim, xlim, '--k')
    #
#     hbins = numpy.linspace(-4, 4, num=161)
#     fig = pyplot.figure()
#     pyplot.subplot(121, aspect='equal')
#     pyplot.hist2d(vl_A, vl_B, hbins, cmap='gray_r')
#     pyplot.subplot(122, aspect='equal')
#     pyplot.hist2d(vc_A, vc_B, hbins, cmap='gray_r')
    #
#     pyplot.show()

    # as dataframe
    data_A = pandas.DataFrame({'vl':vl_A, 'vc':vc_A}, index=t_A)
    data_B = pandas.DataFrame({'vl':vl_B, 'vc':vc_B}, index=t_B)
    # detect outliers?
    data_A = data_A.resample('1S').mean().dropna()
    data_B = data_B.resample('1S').mean().dropna()
    data_A = rolling_median_test(data_A)
    data_B = rolling_median_test(data_B)
    idx = data_A.index.intersection(data_B.index)

    # RMS difference
    vl_rms = numpy.sqrt(numpy.mean((data_A['vl'][idx] - data_B['vl'][idx])**2))
    vc_rms = numpy.sqrt(numpy.mean((data_A['vc'][idx] - data_B['vc'][idx])**2))
    # linear regression
    vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(data_A['vl'][idx], data_B['vl'][idx])
    vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(data_A['vc'][idx], data_B['vc'][idx])
    print 'Filtered (1S) stats:'.format(resample)
    print '\tRMS: {:.2f}, {:.2f}'.format(vl_rms, vc_rms)
    print '\tslopes: {:.2f}, {:.2f}'.format(vl_slope, vc_slope)
    print '\toffsets: {:.2f}, {:.2f}'.format(vl_offset, vc_offset)
    print '\tr^2: {:.2f}, {:.2f}'.format(vl_r**2, vc_r**2)
    print '\tp: {:.2f}, {:.2f}'.format(vl_p, vc_p)
    # plot
    fig = pyplot.figure()
    ax1 = pyplot.subplot(121, aspect='equal')
    pyplot.plot(data_A['vl'][idx], data_B['vl'][idx], '.', alpha=0.5)
    xlim = ax1.get_xlim()
    pyplot.plot(xlim, xlim, '--k')
    ax2 = pyplot.subplot(122, aspect='equal')
    pyplot.plot(data_A['vc'][idx], data_B['vc'][idx], '.', alpha=0.5)
    xlim = ax2.get_xlim()
    pyplot.plot(xlim, xlim, '--k')


    # running mean
    data_Ar = data_A.resample(resample).mean()
    data_Br = data_B.resample(resample).mean()

    # RMS difference
    vl_rms = numpy.sqrt(numpy.mean((data_Ar['vl'] - data_Br['vl'])**2))
    vc_rms = numpy.sqrt(numpy.mean((data_Ar['vc'] - data_Br['vc'])**2))
    # linear regression
    vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(data_Ar['vl'], data_Br['vl'])
    vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(data_Ar['vc'], data_Br['vc'])
    print 'Resampled ({}) stats:'.format(resample)
    print '\tRMS: {:.2f}, {:.2f}'.format(vl_rms, vc_rms)
    print '\tslopes: {:.2f}, {:.2f}'.format(vl_slope, vc_slope)
    print '\toffsets: {:.2f}, {:.2f}'.format(vl_offset, vc_offset)
    print '\tr^2: {:.2f}, {:.2f}'.format(vl_r**2, vc_r**2)
    print '\tp: {:.2f}, {:.2f}'.format(vl_p, vc_p)
    # plot
    fig = pyplot.figure()
    ax1 = pyplot.subplot(121, aspect='equal')
    pyplot.plot(data_Ar['vl'], data_Br['vl'], '.', alpha=0.5)
    xlim = ax1.get_xlim()
    pyplot.plot(xlim, xlim, '--k')
    ax2 = pyplot.subplot(122, aspect='equal')
    pyplot.plot(data_Ar['vc'], data_Br['vc'], '.', alpha=0.5)
    xlim = ax2.get_xlim()
    pyplot.plot(xlim, xlim, '--k')

    pyplot.show()

def compare_ADV(resample='10S'):

    # load ADV
    print 'Loading ADV...'
    t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130413.txt',
                                                rotate=-10.)
    data_ADV = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)

    # load estimates
    print 'Loading estimates...'
    flist_30 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/timeseries_30m/timeseries_30m_20140313_*_probe.json'))
    #flist_60 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/timeseries_60m/timeseries_60m_20140313_*_probe.json'))
    #flist_125 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/swash_125m/20140313/20140313_15/swash_125m_20140313_15*_probe.json'))
    t_estim, vl_estim, vc_estim = load_OF_timeseries(flist_30, rotate=-10.)
    data_estim = pandas.DataFrame({'vl':vl_estim, 'vc':vc_estim}, index=t_estim)

    # crop
#     data_estim = data_estim['2014-03-13 15:00':'2014-03-13 16:00']
#     data_ADV = data_ADV['2014-03-13 15:00':'2014-03-13 16:00']

    # steps:
    # - remove outliers (10%)
    data_ADV = rolling_median_test2(data_ADV)
    data_estim = rolling_median_test2(data_estim)
    # - resample at 1 Hz
    data_ADV = data_ADV.resample('1S').mean()
    data_estim = data_estim.resample('1S').mean()
    # - smooth out over 3 pts
    data_ADVr = data_ADV.rolling(3, win_type='triang').mean()
    data_estimr = data_estim.rolling(3, win_type='triang').mean()
    # - resample
    data_ADVr = data_ADVr.resample(resample).mean()
    data_estimr = data_estimr.resample(resample).mean()

    # filter / average
#     data_ADV = rolling_median_test2(data_ADV)
#     data_ADVr = data_ADV.resample(resample).mean()
#     data_estim = rolling_median_test2(data_estim)
#     data_estimr = data_estim.resample(resample).mean()

    # intersect index of non-null data
    idx = data_estimr.dropna().index.intersection(data_ADVr.dropna().index)
    vl_rms = numpy.sqrt(numpy.mean((data_ADVr['vl'][idx] - data_estimr['vl'][idx])**2))
    vc_rms = numpy.sqrt(numpy.mean((data_ADVr['vc'][idx] - data_estimr['vc'][idx])**2))
    vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(data_ADVr['vl'][idx], data_estimr['vl'][idx])
    vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(data_ADVr['vc'][idx], data_estimr['vc'][idx])
    print '{:.2f} ({:.2f}, {:.2f}, {:.2f}); {:.2f} ({:.2f}, {:.2f}, {:.2f})'.format(
        vl_r, vl_slope, vl_offset, vl_rms, vc_r, vc_slope, vc_offset, vc_rms)


#     print len(idx)
#     u_A = data_estimr['vl'][idx]
#     v_A = data_estimr['vc'][idx]
#     dir_A = numpy.rad2deg(numpy.arctan2(v_A, u_A))
#     x_B = data_ADVr['vl'][idx]#data_estimr['vl'][idx]
#     y_B = data_ADVr['vc'][idx]#data_estimr['vc'][idx]
#     angles = range(-90, 90, 2)
#     r_l = []
#     r_c = []
#     for a in angles:
#         cosa = numpy.cos(numpy.deg2rad(a))
#         sina = numpy.sin(numpy.deg2rad(a))
#         u_B = cosa*x_B - sina*y_B
#         v_B = sina*x_B + cosa*y_B
#         dir_B = numpy.rad2deg(numpy.arctan2(v_B, u_B))
#         vl_rms = numpy.sqrt(numpy.mean((u_A - u_B)**2))
#         vc_rms = numpy.sqrt(numpy.mean((v_A - v_B)**2))
#         vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(u_A, u_B)
#         vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(v_A, v_B)
#         dir_slope, dir_offset, dir_r, dir_p, dir_stderr = stats.linregress(dir_A, dir_B)
#         print '{:03d}: {:.2f} ({:.2f}, {:.2f}, {:.2f}); {:.2f} ({:.2f}, {:.2f}, {:.2f})'.format(
#             a, vl_r, vl_slope, vl_offset, vl_rms, vc_r, vc_slope, vc_offset,
#             vc_rms)
#         r_l.append(vl_r)
#         r_c.append(vc_r)
#     pyplot.figure()
#     pyplot.plot(angles, r_l, '-k')
#     pyplot.plot(angles, r_c, '-r')

    # plot time-series
    fig1 = pyplot.figure()
    ax1 = pyplot.subplot(211, ylabel='long-shore')
    ax2 = pyplot.subplot(212, ylabel='cross-shore')
    #ax3 = pyplot.subplot(313, ylabel='direction')
    for ax, var in zip([ax1, ax2], ['vl', 'vc']):
        for data, color, label in zip([data_ADVr, data_estimr], ['g', 'b'], ['ADV', 'OF']):
            ax.plot(idx, data[var][idx], color=color, label=label)
            ax.axhline(0., ls='--', color='k')
    #for data, color, label in zip([data_ADVr, data_estimr], ['g', 'b'], ['ADV', 'OF']):
    #    dir = numpy.rad2deg(numpy.arctan2(data['vc'][idx], data['vl'][idx]))
    #    ax3.plot(idx, dir, '+', color=color, label=label)
    #ax3.set_ylim(-90., 90.)
    # set nice time formatter
    for ax in [ax1, ax2,]:
        ax.xaxis.set_major_locator(dates.MinuteLocator([0, 30]))
        ax.xaxis.set_minor_locator(dates.MinuteLocator([10, 20, 40, 50]))
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax2.set_xlabel('Video time - {:%Y-%m-%d}'.format(data_estimr.index[0]))
    ax1.legend(ncol=2, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.01))

    # plot scatter
    fig = pyplot.figure()
    ax1 = pyplot.subplot(121, aspect='equal',)
    pyplot.plot(data_ADVr['vl'][idx], data_estimr['vl'][idx], '.', alpha=0.15)
    xlim = ax1.get_xlim()
    pyplot.plot(xlim, xlim, '--k')
    ax2 = pyplot.subplot(122, aspect='equal')
    pyplot.plot(data_ADVr['vc'][idx], data_estimr['vc'][idx], '.', alpha=0.15)
    xlim = ax2.get_xlim()
    pyplot.plot(xlim, xlim, '--k')
    pyplot.show()

def compare_ADV_2(resample='10S', rotate=-10.):
    """ Interpolates ADV data at OF times
    """
    # load ADV
    print 'Loading ADV...'
    t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130413.txt',
                                                rotate=rotate)
    tf_ADV = dates.date2num(t_ADV)

    # load estimates
    print 'Loading v2 estimates...'
    flist_30 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/timeseries_30m/timeseries_30m_20140313_*_probe.json'))
    flist_60 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/timeseries_60m/timeseries_60m_20140313_*_probe.json'))
    #flist_125 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/swash_125m/20140313/20140313_15/swash_125m_20140313_15*_probe.json'))
    t_estim, vl_estim, vc_estim = load_OF_timeseries(flist_60, rotate=rotate)
    tf_estim = dates.date2num(t_estim)

    # make pandas objects
    data_ADVi = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)
    data_estim = pandas.DataFrame({'vl':vl_estim, 'vc':vc_estim}, index=t_estim)

    # filter out
    tau = 0
    if tau>0:
        data_estim = rolling_median_test(data_estim, tau=tau)
        data_ADVi = rolling_median_test(data_ADVi, tau=tau)

    # resample/average
    data_ADVi = data_ADVi.resample(resample).mean().dropna()
    data_estim = data_estim.resample(resample).mean().dropna()

    # crop
    tmin = '2014-03-13 15:00'
    tmax = '2014-03-13 16:00'
    data_estim = data_estim[tmin:tmax]
    data_ADVi = data_ADVi[tmin:tmax]

    #stats
    idx = data_estim.index.intersection(data_ADVi.index)
    print 'Statistics...'
    vl_rms = numpy.sqrt(numpy.mean((data_ADVi['vl'][idx] - data_estim['vl'][idx])**2))
    vc_rms = numpy.sqrt(numpy.mean((data_ADVi['vc'][idx] - data_estim['vc'][idx])**2))
    vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(data_ADVi['vl'][idx], data_estim['vl'][idx])
    vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(data_ADVi['vc'][idx], data_estim['vc'][idx])
    print 'r2={:.2f} ({:.2f}, {:.2f}, {:.2f}); r2={:.2f} ({:.2f}, {:.2f}, {:.2f})'.format(
        vl_r**2, vl_slope, vl_offset, vl_rms, vc_r**2, vc_slope, vc_offset, vc_rms)

    # plot time-series
    fig1 = pyplot.figure()
    ax1 = pyplot.subplot(211, ylabel='long-shore')
    ax2 = pyplot.subplot(212, ylabel='cross-shore')
    #ax3 = pyplot.subplot(313, ylabel='direction')
    for ax, var in zip([ax1, ax2], ['vl', 'vc']):
        for data, color, label, lw in zip([data_ADVi, data_estim],
                                          ['.6', 'b'], ['ADV', 'OF'], [3., 1.5]):
            ax.plot(data.index, data[var], color=color, label=label, lw=lw)
            ax.axhline(0., ls='--', color='k')
    # set nice time formatter
    for ax in [ax1, ax2,]:
        ax.xaxis.set_major_locator(dates.MinuteLocator([0, 30]))
        ax.xaxis.set_minor_locator(dates.MinuteLocator([10, 20, 40, 50]))
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax2.set_xlabel('Video time - {:%Y-%m-%d}{}{}'.format(
        data_estim.index[0],
        (' - {} average'.format(resample) if resample is not None else ''),
        (' - tau={} median test'.format(tau) if tau>0 else ''),
        ))
    ax1.legend(ncol=2, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.01))

    # plot scatter
    fig = pyplot.figure()
    ax1 = pyplot.subplot(121, aspect='equal',)
    pyplot.plot(data_ADVi['vl'][idx], data_estim['vl'][idx], '.', alpha=0.15)
    ax2 = pyplot.subplot(122, aspect='equal')
    pyplot.plot(data_ADVi['vc'][idx], data_estim['vc'][idx], '.', alpha=0.15)
    for ax in [ax1, ax2]:
        xlim = ax.get_xlim()
        ax.plot(xlim, xlim, '-k')
        ax.set_xlim(xlim)
        ax.set_ylim(xlim)
    pyplot.show()

def compare_ADV_3(rotate=-10.):

    print 'Reference rotation =', rotate
    # load ADV
    print 'Loading ADV...'
    t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130413_refEastNorth.txt', rotate=rotate)
    tf_ADV = dates.date2num(t_ADV)

    # load estimates
    print 'Loading v2 estimates...'
    flist_30 = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_30m/timeseries_30m_20140313_*_probe.json')))
    flist_60 = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_60m/timeseries_60m_20140313_*_probe.json')))
    t_estim, vl_estim, vc_estim = load_OF_timeseries(flist_30, rotate=rotate)

    # make pandas objects
    data_ADVi = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)
    data_estim = pandas.DataFrame({'vl':vl_estim, 'vc':vc_estim}, index=t_estim)

    # resample (and average), drop missing data
    resample = '1S'
    data_ADVi = data_ADVi.resample(resample).mean().dropna()
    data_estim = data_estim.resample(resample).mean().dropna()

    # crop
    tmin = '2014-03-13 12:30'
    tmax = '2014-03-13 17:00'
    data_estim = data_estim[tmin:tmax]
    data_ADVi = data_ADVi[tmin:tmax]
    zmin = 1000
    zmax = zmin + 900

    # spectrum
#     win = 300
#     dt_ADVi = (data_ADVi.index[1] - data_ADVi.index[0]).total_seconds()
#     nfft_ADVi = int(win/dt_ADVi)
#     dt_estim = (data_estim.index[1] - data_estim.index[0]).total_seconds()
#     nfft_estim = int(win/dt_estim)
#     print dt_ADVi, nfft_ADVi
#     print dt_estim, nfft_estim
#     f_ADVi, Psd_ADVi = signal.periodogram(data_ADVi['vl'], fs=1./dt_ADVi, nfft=nfft_ADVi)
#     f_estim, Psd_estim = signal.periodogram(data_estim['vl'], fs=1./dt_estim, nfft=nfft_estim)
#
#     # plot spectrum
#     fig = pyplot.figure()
#     ax = fig.add_subplot(111, xlabel='f [Hz]', ylabel='PSD [m^2/s^2]',
#                          xscale='log', yscale='log')
#     ax.plot(f_ADVi[1:], Psd_ADVi[1:])
#     ax.plot(f_estim[1:], Psd_estim[1:])

    # plot scatter

    dpi = 90.
    fig1 = pyplot.figure(figsize=(1280./dpi, 720./dpi))
    ax1 = pyplot.subplot(121, xlabel='ADV V [m/s]', ylabel='video V [m/s]',
                         title='Reference rotation w.r.t. East/North: {} degree'.format(rotate),
                         aspect='equal', xlim=[-2.5, 2.5], ylim=[-2.5, 2.5],
                         )
    ax2 = pyplot.subplot(122, xlabel='ADV U [m/s]', ylabel='video U [m/s]',
                         aspect='equal', xlim=[-2.5, 2.5], ylim=[-2.5, 2.5],
                         )
#     edges = numpy.linspace(-2., 2., 101)
    scatter_resample = '7S'
    tmp1 = data_estim.resample(scatter_resample).mean()
    tmp2 = data_ADVi.resample(scatter_resample).mean()
    in_both = tmp1.index.intersection(tmp2.index)
    print len(in_both)
    for var, ax in zip(['vl', 'vc'], [ax1, ax2]):
        ax.scatter(tmp1[var][in_both], tmp2[var][in_both], 2.,
                   alpha=0.05)
        ax.plot([-2.5, 2.5], [-2.5, 2.5], '--k')
        # note: first var is vertical, second horizontal (cf doc)
#         hist, _, _ = numpy.histogram2d(data_estim[var][in_both], data_ADVi[var][in_both], [edges, edges])
#         print hist.shape, hist.dtype
#         p=ax.imshow(hist, origin='lower', cmap='gray_r')#, vmin=0., vmax=0.1)#,
#                   extent=[edges[0], edges[-1], edges[0], edges[-1]])
#     pyplot.colorbar(p)
    pyplot.show()
#     return
    # plot time-series
    for var, label, ylim in zip(['vl', 'vc'],
                                ['long-shore V', 'cross-shore U'],
                                [[-1., 1.5], [-1., 1.5]]):
        fig1 = pyplot.figure(figsize=(1280./dpi, 720./dpi))
        ax1 = pyplot.subplot(211, ylabel='{} [m/s]'.format(label),
                             title='Reference rotation w.r.t. East/North: {} degree'.format(rotate),
                             )
        ax2 = pyplot.subplot(212, ylabel='{} [m/s]'.format(label),
                             xlabel='Video time - {:%Y-%m-%d} UTC'.format(data_estim.index[0]),
                             )
        for data, color, label, lw, zorder in zip([data_estim, data_ADVi],
                                                  ['r', '.4'], ['OF', 'ADV'],
                                                  [1., 1.], [2, 1]):
            # the enveloppe
            # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html
            vmean = numpy.mean(data[var])
            ht = signal.hilbert(data[var] - vmean)
            data['envp'] = vmean + numpy.abs(ht)
            data['envm'] = vmean - numpy.abs(ht)
            # rolling averages
            tmp1 = data.rolling(121, center=False, min_periods=1).mean() #120 pt (2-min) rolling window
            tmp2 = data.rolling(2, center=False, min_periods=1).mean() #2pt (2 s) rolling window
            # curves
            ax1.plot(tmp1.index, tmp1[var], color=color, label=label, lw=2*lw, zorder=10+zorder)
    #         ax1.plot(tmp1.index, tmp1['env'], color=color, label=label, lw=lw, zorder=10+zorder)
            ax1.fill_between(tmp1.index, tmp1['envm'], tmp1['envp'],
                             facecolor=color, edgecolor=color, alpha=0.25, zorder=zorder)
            ax2.plot(tmp2.index, tmp2[var], color=color, label=label, lw=lw, zorder=zorder)
        ax1.axhline(0., ls='--', color='k')
        ax2.axhline(0., ls='--', color='k')
        # formatters
        ax1.xaxis.set_major_locator(dates.MinuteLocator([0, 30]))
        ax1.xaxis.set_minor_locator(dates.MinuteLocator([10, 20, 40, 50]))
        ax1.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
        ax2.xaxis.set_major_locator(dates.MinuteLocator(range(0,60,2)))
        ax2.xaxis.set_minor_locator(dates.MinuteLocator(range(60)))
        ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
        # limits
        ax1.set_ylim(ylim)
        ax1.axvspan(data_estim.index[zmin], data_estim.index[zmax], zorder=0,
                    facecolor='.8', edgecolor='none')
        ax2.set_xlim(data_estim.index[zmin], data_estim.index[zmax])
        ax2.set_ylim(-2., 2.)
        # adjust
        pyplot.subplots_adjust(left=.06, bottom=.07, right=.98, top=.96)
        outfile = '../figures/misc/timeseries_{}_rot{}{}.png'.format(
            var, ('p' if rotate>0 else 'm'), int(abs(rotate)))
        fig1.savefig(outfile,
                     dpi=dpi)
        print 'saved', outfile
    pyplot.show()

#########

def compare_estim(resample='10S'):
    """Compare OF estimates from v1 (.mat) with v2 estimates (.json)
    """
    # load V2 estimates
    print 'Loading estimates...'
    flist_30 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/timeseries_30m/timeseries_30m_20140313_*_probe.json'))
    flist_60 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/timeseries_60m/timeseries_60m_20140313_*_probe.json'))
    #flist_125 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/swash_125m/20140313/20140313_15/swash_125m_20140313_15*_probe.json'))
    t_estim, vl_estim, vc_estim = load_OF_timeseries(flist_60, rotate=0.)

    # load V1 estimates
    olddata = scio.loadmat('resources/series_median_5_20140313-mix.mat')
    vl_old = numpy.squeeze(olddata['vx'])
    vc_old = numpy.squeeze(olddata['vy'])
    t_old = [datetime.datetime.strptime(ts, '%Y-%m-%d %H:%M:%S.%f') for ts in olddata['t']]

    # load ADV
    print 'Loading ADV...'
    t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130413.txt',
                                                rotate=0.)

    # make pandas objects
    data_old = pandas.DataFrame({'vl':vl_old, 'vc':vc_old}, index=t_old)
    data_estim = pandas.DataFrame({'vl':vl_estim, 'vc':vc_estim}, index=t_estim)
    data_ADV = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)

    # filter
    tau = 0.
    if tau>0:
        data_old = rolling_median_test2(data_old, tau=tau)
        data_estim = rolling_median_test2(data_estim, tau=tau)
        data_ADV = rolling_median_test2(data_ADV, tau=tau)

    # resample/average
    if resample:
        data_old = data_old.resample(resample).mean().dropna()
        data_estim = data_estim.resample(resample).mean().dropna()
        data_ADV = data_ADV.resample(resample).mean().dropna()

    # crop
    tmin = '2014-03-13 15:00'
    tmax = '2014-03-13 16:00'
    data_estim = data_estim[tmin:tmax]
    data_old = data_old[tmin:tmax]
    data_ADV = data_ADV[tmin:tmax]

    #stats
    idxOF = data_estim.index.intersection(data_old.index)
    idxADV = data_ADV.index.intersection(idxOF)
    print 'Statistics, old OF vs new OF...'
    vl_rms = numpy.sqrt(numpy.mean((data_old['vl'][idxOF] - data_estim['vl'][idxOF])**2))
    vc_rms = numpy.sqrt(numpy.mean((data_old['vc'][idxOF] - data_estim['vc'][idxOF])**2))
    vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(data_old['vl'][idxOF], data_estim['vl'][idxOF])
    vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(data_old['vc'][idxOF], data_estim['vc'][idxOF])
    print 'r2={:.2f} ({:.2f}, {:.2f}, {:.2f}); r2={:.2f} ({:.2f}, {:.2f}, {:.2f})'.format(
        vl_r**2, vl_slope, vl_offset, vl_rms, vc_r**2, vc_slope, vc_offset, vc_rms)

    # plot time-series
    fig1 = pyplot.figure()
    ax1 = pyplot.subplot(211, ylabel='eastern')
    ax2 = pyplot.subplot(212, ylabel='northern')
    #ax3 = pyplot.subplot(313, ylabel='direction')
    for ax, var in zip([ax1, ax2], ['vl', 'vc']):
        for data, idx, color, label, zorder, lw, ls in zip(
            [data_old, data_estim, data_ADV],
            [idxOF, idxOF, idxADV],
            ['g', 'c', '0.8'],
            ['old OF', 'new OF', 'ADV'],
            [3,4,2], [1.5, 1.5, 3.], ['-', '-', '-']):
            ax.plot(idx, data[var][idx], color=color, label=label, zorder=zorder,
                    lw=lw, ls=ls)
            ax.axhline(0., ls='--', color='k')
    # set nice time formatter
    for ax in [ax1, ax2,]:
        ax.xaxis.set_major_locator(dates.MinuteLocator([0, 30]))
        ax.xaxis.set_minor_locator(dates.MinuteLocator([10, 20, 40, 50]))
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax2.set_xlabel('Video time - {:%Y-%m-%d}{}{}'.format(
        data_estim.index[0],
        (' - {} average'.format(resample) if resample is not None else ''),
        (' - tau={} median test'.format(tau) if tau>0 else ''),
        ))
    ax1.legend(ncol=3, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.01))

    # plot scatter
    fig = pyplot.figure()
    ax1 = pyplot.subplot(121, aspect='equal',)
    pyplot.plot(data_old['vl'][idxOF], data_estim['vl'][idxOF], '.', alpha=0.15)
    xlim = ax1.get_xlim()
    pyplot.plot(xlim, xlim, '--k')
    ax2 = pyplot.subplot(122, aspect='equal')
    pyplot.plot(data_old['vc'][idxOF], data_estim['vc'][idxOF], '.', alpha=0.15)
    xlim = ax2.get_xlim()
    pyplot.plot(xlim, xlim, '--k')
    pyplot.show()

def compare_old(resample=None, tau=0):

    # load ADV
    print 'Loading ADV...'
    imin = 1
    imax = 460800 # the last index of 2014-03-13 data
    rawADV = scio.loadmat('resources/vel_13-16mars2014_XYZ.mat')
    t0_ADV = datetime.datetime(2014, 3, 13, 8, 0, 0) # excluding the 189s shift (3 min 9 s)
    t_ADV = [t0_ADV+datetime.timedelta(seconds=s) for s in rawADV['A'][imin:imax+1,0]]
    UADV = 0.01*rawADV['A'][imin:imax+1,1]
    VADV = 0.01*rawADV['A'][imin:imax+1:,2]
    tf_ADV = dates.date2num(t_ADV)

    # load V1 estimates
    print 'Loading v1 estimates...'
    estimdata = scio.loadmat('resources/series_median_5_20140313-mix.mat')
    vx_estim = numpy.squeeze(estimdata['vx'])
    vy_estim = numpy.squeeze(estimdata['vy'])
    t_estim = [datetime.datetime.strptime(ts, '%Y-%m-%d %H:%M:%S.%f') for ts in estimdata['t']]
    tf_estim = dates.date2num(t_estim)

    # interpolate ADV at OF times
    dt = (0.5*377)/(24.*3600.)
    print 'Interpolation...'
    interpolator = interpolate.InterpolatedUnivariateSpline(tf_ADV, UADV)
    UADVi = interpolator(tf_estim - dt)
    interpolator = interpolate.InterpolatedUnivariateSpline(tf_ADV, VADV)
    VADVi = interpolator(tf_estim - dt)
    t_ADVi = t_estim

    # rotate the video estimates
    alpha = 25.;
    cos_alpha = numpy.cos(numpy.deg2rad(alpha))
    sin_alpha = numpy.sin(numpy.deg2rad(alpha))
    # NOTE: this rotation is WRONG!!!!!!
    UVID = -vy_estim*cos_alpha + vx_estim*sin_alpha;
    VVID = vx_estim*cos_alpha - vy_estim*sin_alpha;
    speed_estim = numpy.sqrt(vx_estim**2 + vy_estim**2)
    speedVID = numpy.sqrt(UVID**2 + VVID**2)
    print speed_estim/speedVID

    # make pandas objects
    data_ADVi = pandas.DataFrame({'u':UADVi, 'v':VADVi}, index=t_ADVi)
    data_estim = pandas.DataFrame({'u':UVID, 'v':VVID}, index=t_estim)

    # resample/average
    data_ADVi = data_ADVi.rolling(2).mean()
    data_estim = data_estim.rolling(2).mean()
    data_ADVi = data_ADVi.rolling(5).mean()
    data_estim = data_estim.rolling(5).mean()
    #data_ADVi = data_ADVi.resample(resample).mean().dropna()
    #data_estim = data_estim.resample(resample).mean().dropna()

    # crop
    tmin = '2014-03-13 15:00'
    tmax = '2014-03-13 16:00'
    data_estim = data_estim[tmin:tmax]
    data_ADVi = data_ADVi[tmin:tmax]

    # load the pseudo "scatterplot" data
#     scatterdata = scio.loadmat('resources/Data_Scatter_PLot_ADV_OF_GPP2.mat');
#     tmp2 = numpy.squeeze(scatterdata['V_OF'])
#     tmp1 = numpy.squeeze(data_estim.as_matrix(['v',]))
#     print tmp1.shape, tmp1
#     print tmp2.shape, tmp2
#     x1 = numpy.linspace(0, 1., tmp1.size)
#     x2 = numpy.linspace(0, 1., tmp2.size)
#     pyplot.figure()
#     pyplot.plot(x1[:1000], tmp1[:1000], 'k')
#     pyplot.plot(x2[:1000], tmp2[:1000], 'r')
#     pyplot.show()
#     return

    #stats
    idx = data_estim.index.intersection(data_ADVi.index)

    tmp1 = data_ADVi['u'][idx].values
    tmp2 = data_estim['u'][idx].values
    print tmp1.shape, tmp1
    print tmp2.shape, tmp2
    vl_rms =  numpy.sqrt(numpy.mean((tmp1 - tmp2)**2))
    vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(tmp1, tmp2)
    print 'r2={:.2f} ({:.2f}, {:.2f}, {:.2f});'.format(
        vl_r**2, vl_slope, vl_offset, vl_rms)

    print 'Statistics...'
    vl_rms = numpy.sqrt(numpy.mean((data_ADVi['u'][idx] - data_estim['u'][idx])**2))
    vc_rms = numpy.sqrt(numpy.mean((data_ADVi['v'][idx] - data_estim['v'][idx])**2))
    vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(data_ADVi['u'][idx], data_estim['u'][idx])
    vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(data_ADVi['v'][idx], data_estim['v'][idx])
    print 'r2={:.2f} ({:.2f}, {:.2f}, {:.2f}); r2={:.2f} ({:.2f}, {:.2f}, {:.2f})'.format(
        vl_r**2, vl_slope, vl_offset, vl_rms, vc_r**2, vc_slope, vc_offset, vc_rms)

    # plot time-series
    fig1 = pyplot.figure()
    ax1 = pyplot.subplot(211, ylabel='u')
    ax2 = pyplot.subplot(212, ylabel='v')
    for ax, var in zip([ax1, ax2], ['u', 'v']):
        for data, color, label, lw in zip([data_ADVi, data_estim],
                                          ['.6', 'b'], ['ADV', 'OF'], [3., 1.]):
            ax.plot(data.index, data[var], color=color, label=label, lw=lw)
            ax.axhline(0., ls='--', color='k')
    # set nice time formatter
    for ax in [ax1, ax2,]:
        ax.xaxis.set_major_locator(dates.MinuteLocator([0, 30]))
        ax.xaxis.set_minor_locator(dates.MinuteLocator([10, 20, 40, 50]))
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax2.set_xlabel('Video time - {:%Y-%m-%d}{}{}'.format(
        data_estim.index[0],
        (' - {} average'.format(resample) if resample is not None else ''),
        (' - tau={} median test'.format(tau) if tau>0 else ''),
        ))
    ax1.legend(ncol=2, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.01))

    # plot scatter
    fig = pyplot.figure()
    ax1 = pyplot.subplot(121, aspect='equal',)
    pyplot.plot(data_ADVi['u'][idx], data_estim['v'][idx], '.', alpha=0.15)
    ax2 = pyplot.subplot(122, aspect='equal')
    pyplot.plot(data_ADVi['v'][idx], data_estim['v'][idx], '.', alpha=0.15)
    for ax in [ax1, ax2]:
        xlim = ax.get_xlim()
        ax.plot(xlim, xlim, '-k')
        ax.set_xlim(xlim)
        ax.set_ylim(xlim)
    pyplot.show()

def export_series(outdir=None):
    """Export the time series from json to mat files.
    """
    # the files to load
    flist_30 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/timeseries_30m/timeseries_30m_20140313_*_probe.json'))
    flist_60 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/timeseries_60m/timeseries_60m_20140313_*_probe.json'))

    # for each set
    for fl, pr in zip([flist_30, flist_60],
                      [PARAMS_COMP30, PARAMS_COMP60]):
        # load estim (NO ROTATION!)
        t_estim, vx_estim, vy_estim = load_OF_timeseries(flist_60, rotate=0.)
        # time as str
        tstr_estim = [t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in t_estim]
        # the data to be archived
        data = {'author': 'Pierre DERIAN',
                'created_by': '{} on {}'.format(__file__, datetime.datetime.now()),
                'description': 'Estimated velocity field by Typhoon from images gridded at resolution {} m/pixel and low-pass median_filtered with a {}-m filter. Measures are given in a reference rotated of {} degree CCW w.r.t. the Eastern / Northern reference. Values were averaged in a {} m radius circle center at (x,y)=({},{}) UTM.'.format(pr['resolution'],
             pr['median_length'],
             pr['rotation'],
             AVG_PROBE['r'],
             AVG_PROBE['x'],
             AVG_PROBE['y'],),
                'label': pr['label'],
                'vx': vx_estim,
                'vx_descr': 'velocities along x [m/s]',
                'vy': vy_estim,
                'vy_descr': 'velocities along y [m/s]',
                't': tstr_estim,
                't_descr': 'velocities timestamp, format "%Y-%m-%d %H:%M:%S.%f"',
                }
        if outdir is not None:
            outfile = os.path.join(outdir, '{}.mat'.format(pr['label']))
            scio.savemat(outfile, data)
            print 'saved', outfile

def test_timeshift(resample='1S'):

    # load estimates
    print 'Loading estimates...'
    flist_30 = sorted(glob.glob('../data/plage/estim/timeseries_30m/timeseries_30m_20140313_15*_probe.json'))
    flist_60 = sorted(glob.glob('../data/plage/estim/timeseries_60m/timeseries_60m_20140313_15*_probe.json'))
    flist_125 = sorted(glob.glob('/Volumes/LaCie_Mac/pderian/data_GPP/Estim/swash_125m/20140313/20140313_15/swash_125m_20140313_15*_probe.json'))
    t_estim, vl_estim, vc_estim = load_timeseries(flist_30, rotate=0)
    data_estim = pandas.DataFrame({'vl':vl_estim, 'vc':vc_estim}, index=t_estim)
    print '\trolling mean...'
    data_estimr = data_estim.resample(resample).mean()

    # load ADV
    print 'Loading ADV...'
    datestr_ADV, timestr_ADV, u_ADV, v_ADV, w_ADV = numpy.loadtxt(
        'resources/ADVdata_20130413_NOSHIFT.txt',
        dtype=numpy.dtype('a10, a13, f8, f8, f8'),
        unpack=True,
        )
    print '\tgenerating times...'
    t_ADV = [datetime.datetime.strptime(d+' '+t, '%Y-%m-%d %H:%M:%S.%f')
             for d, t in zip(datestr_ADV, timestr_ADV)]

    shifts = numpy.linspace(185., 195., 21)
    r = []
    for s in shifts:
        dt = datetime.timedelta(seconds=s)
        tmp_t = [t+dt for t in t_ADV]
        data_ADV = pandas.DataFrame({'vl':u_ADV, 'vc':-v_ADV}, index=tmp_t)
        data_ADVr = data_ADV.resample(resample).mean()
        idx = data_estimr.index.intersection(data_ADVr.index)
        vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(
            data_estimr['vl'][idx], data_ADVr['vl'][idx])
        print '{:.2f}: {:.2f}'.format(s, vl_r)
        r.append(vl_r)
    pyplot.plot(shifts, r)
    pyplot.show()

def test_misalignment():

    # relevant time window
    tmin = '2014-03-13 12:30'
    tmax = '2014-03-13 17:00'
    tmin_dt = datetime.datetime.strptime(tmin, '%Y-%m-%d %H:%M')
    tmax_dt = datetime.datetime.strptime(tmax, '%Y-%m-%d %H:%M')
    resample_short = '2S' #2 s
    resample_long = '2T' #2 min

    # load ADV, rotate by +25 to restore original ADV reference
    print 'Loading ADV...'
    t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130413_refEastNorth.txt', rotate=25.)
    data_ADV = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)
    data_ADV = data_ADV[tmin:tmax]
    data_ADV = data_ADV.resample(resample_short).mean().dropna()

    # load estimates, rotate by -10 (East/North => long/cross reference)
    print 'Loading v2 estimates...'
    flist_30 = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_30m/timeseries_30m_20140313_*_probe.json')))
    flist_60 = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_60m/timeseries_60m_20140313_*_probe.json')))
    t_estim, vl_estim, vc_estim = load_OF_timeseries(flist_60, rotate=-10.)
    data_estim = pandas.DataFrame({'vl':vl_estim, 'vc':vc_estim}, index=t_estim)
    data_estim = data_estim[tmin:tmax] #crop to given time window
    data_estim = data_estim.resample(resample_short).mean().dropna() #resample at given frequency

    # intersection
    inter_short = data_estim.index.intersection(data_ADV.index)

    # now for each angle
    angles = numpy.arange(-180, 180)
    correl_vl_short = []
    correl_vc_short = []
    correl_vl_long = []
    correl_vc_long = []
    for a in angles:
        # rotate the ADV data
        cosa = numpy.cos(numpy.deg2rad(a))
        sina = numpy.sin(numpy.deg2rad(a))
        data_ADV['vlr'] = cosa*data_ADV['vl'] - sina*data_ADV['vc']
        data_ADV['vcr'] = sina*data_ADV['vl'] + cosa*data_ADV['vc']
        # correlations
        vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(data_ADV['vlr'][inter_short], data_estim['vl'][inter_short])
        vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(data_ADV['vcr'][inter_short], data_estim['vc'][inter_short])
        correl_vl_short.append(vl_r)
        correl_vc_short.append(vc_r)
        # now resample over longer times
        data_ADV_long = data_ADV.resample(resample_long).mean().dropna()
        data_estim_long = data_estim.resample(resample_long).mean().dropna()
        inter_long = data_estim_long.index.intersection(data_ADV_long.index)
        # correlations
        vl_slope, vl_offset, vl_r, vl_p, vl_stderr = stats.linregress(data_ADV_long['vlr'][inter_long], data_estim_long['vl'][inter_long])
        vc_slope, vc_offset, vc_r, vc_p, vc_stderr = stats.linregress(data_ADV_long['vcr'][inter_long], data_estim_long['vc'][inter_long])
        correl_vl_long.append(vl_r)
        correl_vc_long.append(vc_r)

    pyplot.plot(angles, correl_vl_short, '-b', label='long- @{}'.format(resample_short))
    pyplot.plot(angles, correl_vc_short, '-g', label='cross- @{}'.format(resample_short))
    pyplot.plot(angles, correl_vl_long, '--b', label='long- @{}'.format(resample_long))
    pyplot.plot(angles, correl_vc_long, '--g', label='cross- @{}'.format(resample_long))
    pyplot.legend()
    pyplot.show()

if __name__=="__main__":
    #view('../tmp/test_v2/60m_jpg/20140613_0630_info.json')
    test_misalignment()
    #compare_timeseries()
    #compare_ADV_2()
#     compare_ADV_3()
    #test_timeshift()
    #compare_estim()
    #compare_old()
    #export_series('../data')
