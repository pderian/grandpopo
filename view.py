"""Miscellaneous experiments and tests for Grand Popo results post-processing

Left here for posterity.

Written by P. DERIAN 2016-2017
www.pierrederian.net
"""
###
import argparse
import datetime
import glob
import itertools
import json
import os
import sys
###
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
import matplotlib.dates as dates
import numpy
import pandas
import scipy.interpolate as interpolate
import scipy.io as scio
import scipy.signal as signal
import scipy.stats as stats
import skimage.io as skio
import skimage.transform as sktransform
###
from config import *
sys.path.append('/Users/pderian/Documents/Python/Tools')
import inr

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
    median['res'] = numpy.sqrt((data['ux']-median['ux'])**2 + (data['uy']-median['uy'])**2)
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
    resample = '2S'
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
    nwin = 300
    dt = 2.
    _, Psd_vl_ADVi = signal.welch(data_ADVi['vl'], fs=1./dt, nperseg=nwin)
    _, Psd_vl_estim = signal.welch(data_estim['vl'], fs=1./dt, nperseg=nwin)
    f_ADVi, Psd_vc_ADVi = signal.welch(data_ADVi['vc'], fs=1./dt, nperseg=nwin)
    f_estim, Psd_vc_estim = signal.welch(data_estim['vc'], fs=1./dt, nperseg=nwin)
    # plot spectrum
    fig = pyplot.figure()
    ax = fig.add_subplot(111, xlabel='f [Hz]', ylabel='PSD [m^2/s^2]',
                         xscale='log', yscale='log')
    ax.plot(f_ADVi[1:], Psd_vl_ADVi[1:], label='ADV long')
    ax.plot(f_estim[1:], Psd_vl_estim[1:], label='estim long')
    ax.plot(f_ADVi[1:], Psd_vc_ADVi[1:], label='ADV cross')
    ax.plot(f_estim[1:], Psd_vc_estim[1:], label='estim cross')
    pyplot.legend()
    pyplot.show()
    return

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
    pyplot.legend()
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

def ADVspectrum(rotate=-10.):

    print 'Reference rotation =', rotate
    # load ADV
    print 'Loading ADV...'
    t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130413_refEastNorth.txt', rotate=rotate)
    tf_ADV = dates.date2num(t_ADV)
    data_ADV = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)
    # crop
    tmin = '2014-03-13 12:30'
    tmax = '2014-03-13 17:00'
    data_ADV = data_ADV[tmin:tmax]

    ### plot spectrum
    # compute welch's estimates
    dt = 1./8. #same as resample_short
    nwin = int(90*60*(1./dt)) # T min * 60 s/min * dt sample/s
    _, Psd_vl_ADV = signal.welch(data_ADV['vl'], fs=1./dt, nperseg=nwin)
    f_ADV, Psd_vc_ADV = signal.welch(data_ADV['vc'], fs=1./dt, nperseg=nwin)
    # print stuff
    imaxADV = numpy.argmax(Psd_vc_ADV)
    print 'ADV peak: f={:.3f}, period={:.3f}'.format(f_ADV[imaxADV], 1./f_ADV[imaxADV])
    # create figure
    fig0 = pyplot.figure()
    ax0 = fig0.add_subplot(111, xlabel='f (Hz)', ylabel=r'PSD estimate (m2/s)',
                          xscale='log', yscale='log')
    ax0.plot(f_ADV[1:], Psd_vl_ADV[1:], '-', color='.4')
    ax0.plot(f_ADV[1:], Psd_vc_ADV[1:], '-', color='.4', lw=1.5)
    # tune
    ax0.set_xlim(1./900, 8)
    ax0.set_ylim(1e-1, 1.5e2)

    ### test resampling
    data_ADVd = data_ADV.resample('2T').mean()
    data_ADVr = data_ADV.rolling('2T').mean().resample('2T').last()
    # the 2 above are completely equivalent
    fig1 = pyplot.figure()
    pyplot.plot(data_ADVd.index, data_ADVd['vc'], 'k')
    pyplot.plot(data_ADVr.index, data_ADVr['vc'], 'g')
    pyplot.show()

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
        tmp_df = tmp_f*(tmp_x.mean()/phase_speed)
        result.append((tmp_f, tmp_df, tmp_Pxx))
    # regroup the various spectra
    fref = numpy.arange(0., fs, fs/nperseg)
    P0 = numpy.zeros_like(fref);
    nS = numpy.zeros_like(P0);
    for f, df, Pxx in result:
        fds = f-df #doppler-shifted freqs
        ifds = numpy.digitize(fds, fref)
        for j,k in enumerate(ifds):
            # if the destination bin (k) is valid
            if k>0 and k<P0.size:
                P0[k] += Pxx[j]
                nS[k] += 1.
    has_sample = nS>0
    fref = fref[has_sample]
    P0 = P0[has_sample]/nS[has_sample]
    return fref, P0, result

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

def test_cameradt():

    flist_60 = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_60m/timeseries_60m_20140313_*_probe.json')))
    dt = []
    for f in flist_60:
        with open(f) as fp:
            jsondata = json.load(fp)
            dt += jsondata['dt']
    dt = numpy.array(dt)
    print 'camera dt: mean = {:.3f} s, std = {:.3f} s, max = {:.3f}'.format(numpy.mean(dt), numpy.std(dt), dt.max())
    bins = numpy.linspace(0., 1., 101)
    pyplot.hist(dt, bins)
    pyplot.show()

def test_median_vs_raw(rotate=-10.):
    flist_raw = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_30m_raw/timeseries_30m_20140313_130[0-7]_probe.json')))
    flist_med = sorted(glob.glob(os.path.join(ROOT_ESTIM_DIR, 'timeseries_60m/timeseries_60m_20140313_130[0-7]_probe.json')))
    t_raw, vl_raw, vc_raw = load_OF_timeseries(flist_raw, rotate=rotate)
    t_med, vl_med, vc_med = load_OF_timeseries(flist_med, rotate=rotate)
    data_raw = pandas.DataFrame({'vl':vl_raw, 'vc':vc_raw}, index=t_raw)
    data_med = pandas.DataFrame({'vl':vl_med, 'vc':vc_med}, index=t_med)
    data_raw = data_raw.resample('2S').mean().dropna()
    data_med = data_med.resample('2S').mean().dropna()

    pyplot.figure()
    pyplot.plot(data_med['vl'], data_raw['vl'], '+r')
    pyplot.plot(data_med['vc'], data_raw['vc'], '+k')
    pyplot.figure()
    pyplot.subplot(211)
    pyplot.plot(data_med.index, data_med['vl'], 'r')
    pyplot.plot(data_raw.index, data_raw['vl'], 'k')
    pyplot.subplot(212)
    pyplot.plot(data_med.index, data_med['vc'], 'r')
    pyplot.plot(data_raw.index, data_raw['vc'], 'k')
    pyplot.show()

def test_Dopplershift(rotate=-10.):
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

    # resample (and average), drop missing data
    resample_short = '2S' #2 s
    resample_large = '2T' #2 min
    data_ADV = data_ADV.resample(resample_short).mean().dropna()
    data_estim = data_estim.resample(resample_short).mean().dropna()
    # crop
    tmin = '2014-03-13 12:29'
    tmax = '2014-03-13 17:01'
    data_estim = data_estim[tmin:tmax]
    data_ADV = data_ADV[tmin:tmax]

    ### plot spectrum
    # compute welch's estimates
    dt = 2. #must be the same as resample_short
    fs = 1./dt #sampling frequency
    nwin = int(15*60*fs) #15 min x 60s/min x fs sample/s
    f_estim, P_vc_estim, result_vc_estim = welch_dopplershift(data_estim['vc'], fs, nwin)
    f_ADV, P_vc_ADV, result_vc_ADV = welch_dopplershift(data_ADV['vc'], fs, nwin)

    fig = pyplot.figure()
    pyplot.loglog(f_ADV[1:], P_vc_ADV[1:])
    pyplot.loglog(f_estim[1:], P_vc_estim[1:])
    pyplot.show()

##########

def fulltest(rotate=-10.):
    print 'Reference rotation =', rotate

    # no ADV data???!
    # Note: need reference data for 2014-03-12
    if 0:
        # load ADV
        print 'Loading ADV...'
        t_ADV, vl_ADV, vc_ADV = load_ADV_timeseries('resources/ADVdata_20130412_refEastNorth.txt', rotate=rotate)
        tf_ADV = dates.date2num(t_ADV)
        data_ADV = pandas.DataFrame({'vl':vl_ADV, 'vc':vc_ADV}, index=t_ADV)
        #data_ADV = data_ADV.resample('1S').mean().dropna()
        # crop
        tmin = '2014-03-12 11:15' #11h or 12h?
        tmax = '2014-03-12 11:20'
        data_ADV = data_ADV[tmin:tmax]

        pyplot.subplot(211)
        pyplot.plot(data_ADV.index, data_ADV['vl'])
        pyplot.subplot(212)
        pyplot.plot(data_ADV.index, data_ADV['vc'])
        pyplot.show()

    # drone
    ux_drone = []
    uy_drone = []
    t_drone = []
    x_drone = 750 #approximate ADV coords in [px] in rectified half-res uav images
    y_drone = 300
    XX, YY = numpy.meshgrid(range(1144), range(500))
    is_probed = (XX-x_drone)**2 + (YY-y_drone)**2 < 10**2 #average radius of 10 pixel
    imin_drone = 1
    drone_dir = '/Users/pderian/Documents/Data/GrandPopo/data/drone/halfres/estim_median'
    estimlist = [os.path.join(drone_dir, '{:03d}_UV.inr'.format(k)) for k in xrange(imin_drone, 290)]
    masklist = [os.path.join(drone_dir, '{:03d}_Mask.png'.format(k)) for k in xrange(imin_drone, 290)]
    t0 = datetime.datetime(2014, 3, 12, 11, 17, 14, 650000) #same date as frame #98 of shore video 2014-03-12 11h - 16.mp4
    for idx, efile, mfile in zip(range(len(estimlist)), estimlist, masklist):
        ux, uy = inr.readMotion(efile)
        mask = pyplot.imread(mfile)
        # apply to field
        ux = ux*mask
        uy = uy*mask
#         ux_drone.append(-ux[y_drone, x_drone]) #reverse orientation to match that of shore
#         uy_drone.append(uy[y_drone, x_drone])
        ux_drone.append(-ux[is_probed].mean()) #reverse orientation to match that of shore
        uy_drone.append(uy[is_probed].mean())

        t_drone.append(t0 + datetime.timedelta(seconds=idx))
    tf_drone = dates.date2num(t_drone)

    # shorecam
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
        data.append(tmp)
        nfields += tmp['dt'].size
        # generate what would be the input image name
        basename, _ = os.path.splitext(os.path.basename(str(tmp['file'])))
        framenames += ['{}_{:03d}.png'.format(basename, i) for i in xrange(tmp['dt'].size)]
    # "constant" var
    dx = tmp['dx']
    x = tmp['x']
    y = tmp['y']
    XX, YY = numpy.meshgrid(x, y)
    # group fields and times (assuming same shape)
    ux = numpy.vstack(tuple([d['ux'] for d in data]))
    uy = numpy.vstack(tuple([d['uy'] for d in data]))
    dt = numpy.hstack(tuple([d['dt'] for d in data]))
    t_shore = list(itertools.chain.from_iterable([d['t'][:-1].tolist() for d in data]))
    # convert t to datetimes
    t_shore = [datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S.%f') for s in t_shore]
    # conversion factor
    dx_over_dt = dx/dt
    # now probe
    is_probed = (XX-AVG_PROBE['x'])**2 + (YY-AVG_PROBE['y'])**2 < 1.
    ux_shore = []
    uy_shore = []
    for i in xrange(len(dt)):
        tmp_ux = ux[i,:,:]*dx_over_dt[i] # displacements are now in [m/s]
        tmp_uy = uy[i,:,:]*dx_over_dt[i]
        ux_shore.append(tmp_ux[is_probed].mean())
        uy_shore.append(tmp_uy[is_probed].mean())
    tf_shore = dates.date2num(t_shore)

    # as Pandas data
    data_drone = pandas.DataFrame({'ux':ux_drone, 'uy':uy_drone}, index=t_drone)
    data_shore = pandas.DataFrame({'ux':ux_shore, 'uy':uy_shore}, index=t_shore)
    # crop to same time window
    tmin = max(data_drone.index[0], data_shore.index[0])
    tmax = min(data_drone.index[-1], data_shore.index[-1])
    print tmin, tmax
    data_drone = data_drone[tmin:tmax]
    data_shore = data_shore[tmin:tmax]
    # remove outliers
    data_drone = rolling_median_test(data_drone, win='1T', tau=2)
    data_shore = rolling_median_test(data_shore, win='1T', tau=2)
    # rotate
    if rotate:
        cosp = numpy.cos(numpy.deg2rad(rotate))
        sinp = numpy.sin(numpy.deg2rad(rotate))
        data_drone['vl'] = cosp*data_drone['ux'] - sinp*data_drone['uy']
        data_drone['vc'] = sinp*data_drone['ux'] + cosp*data_drone['uy']
        data_shore['vl'] = cosp*data_shore['ux'] - sinp*data_shore['uy']
        data_shore['vc'] = sinp*data_shore['ux'] + cosp*data_shore['uy']

    # rollign mean
    #data_drone = data_drone.rolling('2S').mean()
    #data_shore = data_shore.rolling('2S').mean()


    ### plot data
    ax0 = pyplot.subplot(211)
    pyplot.plot(data_drone.index, data_drone['ux'], 'r')
    pyplot.ylim(-8., 8.)
    ax0b = ax0.twinx()
    pyplot.plot(data_shore.index, data_shore['ux'], 'k')
    pyplot.ylim(-1.5, 1.5)
    ax0.axhline(0., color='.2', ls='--', zorder=0)
    #
    ax1 = pyplot.subplot(212)
    pyplot.plot(data_drone.index, data_drone['uy'], 'r')
    pyplot.ylim(-8., 8.)
    ax1b = ax1.twinx()
    pyplot.plot(data_shore.index, data_shore['uy'], 'k')
    pyplot.ylim(-2., 2.)
    ax1.axhline(0., color='.2', ls='--', zorder=0)

    ### beach reference data
    if rotate:
        pyplot.figure()
        ax0 = pyplot.subplot(211, ylabel='UAV cross (px)')
        pyplot.plot(data_drone.index, data_drone['vc'], 'r')
        pyplot.ylim(-6., 6.)
        ax0b = ax0.twinx()
        ax0b.set_ylabel('Shore cross (m/s)')
        pyplot.plot(data_shore.index, data_shore['vc'], 'k')
        pyplot.ylim(-1.5, 1.5)
        ax0.axhline(0., color='.2', ls='--', zorder=0)
        ax0.spines['left'].set_color('r')
        ax0b.spines['left'].set_color('r')
        ax0.yaxis.label.set_color('r')
        ax0.tick_params(axis='y', which='both', color='r', labelcolor='r')
        #
        ax1 = pyplot.subplot(212, ylabel='UAV along (px)')
        pyplot.plot(data_drone.index, data_drone['vl'], 'r')
        pyplot.ylim(-6., 6.)
        ax1b = ax1.twinx()
        ax1b.set_ylabel('Shore along (m/s)')
        pyplot.plot(data_shore.index, data_shore['vl'], 'k')
        pyplot.ylim(-1.5, 1.5)
        ax1.axhline(0., color='.2', ls='--', zorder=0)
        ax1.spines['left'].set_color('r')
        ax1b.spines['left'].set_color('r')
        ax1.yaxis.label.set_color('r')
        ax1.tick_params(axis='y', which='both', color='r', labelcolor='r')
    ### NORMALIZED data
    pyplot.figure()
    ax0 = pyplot.subplot(211, ylabel='Eastward (normalized)')
    pyplot.plot(data_drone.index, (data_drone['ux']-numpy.mean(data_drone['ux']))/numpy.std(data_drone['ux']), 'r')
    pyplot.plot(data_shore.index, (data_shore['ux']-numpy.mean(data_shore['ux']))/numpy.std(data_shore['ux']), 'k')
    pyplot.ylim(-3., 3.)
    ax0.axhline(0., color='.2', ls='--', zorder=0)
    #
    ax1 = pyplot.subplot(212, ylabel='Northward (normalized)')
    pyplot.plot(data_drone.index, (data_drone['uy']-numpy.mean(data_drone['uy']))/numpy.std(data_drone['uy']), 'r')
    pyplot.plot(data_shore.index, (data_shore['uy']-numpy.mean(data_shore['uy']))/numpy.std(data_shore['uy']), 'k')
    pyplot.ylim(-3., 3.)
    ax1.axhline(0., color='.2', ls='--', zorder=0)

    ### NORMALIZED, beach frame data
    if rotate:
        pyplot.figure()
        ax0 = pyplot.subplot(211, ylabel='Cross (normalized)')
        pyplot.plot(data_drone.index, (data_drone['vl']-numpy.mean(data_drone['vl']))/numpy.std(data_drone['vl']), 'r')
        pyplot.plot(data_shore.index, (data_shore['vl']-numpy.mean(data_shore['vl']))/numpy.std(data_shore['vl']), 'k')
        pyplot.ylim(-3., 3.)
        ax0.axhline(0., color='.2', ls='--', zorder=0)
        #
        ax1 = pyplot.subplot(212, ylabel='Along (normalized)')
        pyplot.plot(data_drone.index, (data_drone['vc']-numpy.mean(data_drone['vc']))/numpy.std(data_drone['vc']), 'r')
        pyplot.plot(data_shore.index, (data_shore['vc']-numpy.mean(data_shore['vc']))/numpy.std(data_shore['vc']), 'k')
        pyplot.ylim(-3., 3.)
        ax1.axhline(0., color='.2', ls='--', zorder=0)


    pyplot.show()

def test_probe(rotate=-10.):

    # shorecam
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
        data.append(tmp)
        nfields += tmp['dt'].size
        # generate what would be the input image name
        basename, _ = os.path.splitext(os.path.basename(str(tmp['file'])))
        framenames += ['{}_{:03d}.png'.format(basename, i) for i in xrange(tmp['dt'].size)]
    # "constant" var
    dx = tmp['dx']
    x = tmp['x']
    y = tmp['y']
    XX, YY = numpy.meshgrid(x, y)
    # group fields and times (assuming same shape)
    ux = numpy.vstack(tuple([d['ux'] for d in data]))
    uy = numpy.vstack(tuple([d['uy'] for d in data]))
    dt = numpy.hstack(tuple([d['dt'] for d in data]))
    t_shore = list(itertools.chain.from_iterable([d['t'][:-1].tolist() for d in data]))
    # convert t to datetimes
    t_shore = [datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S.%f') for s in t_shore]
    # conversion factor
    dx_over_dt = dx/dt
    # now probe
    iProbe = numpy.argmin((y-AVG_PROBE['y'])**2)
    jProbe = numpy.argmin((x-AVG_PROBE['x'])**2)
    is_avg = (XX-AVG_PROBE['x'])**2 + (YY-AVG_PROBE['y'])**2 < 2.
    ux_probe = []
    uy_probe = []
    ux_avg = []
    uy_avg = []
    for i in xrange(len(dt)):
        tmp_ux = ux[i,:,:]*dx_over_dt[i] # displacements are now in [m/s]
        tmp_uy = uy[i,:,:]*dx_over_dt[i]
        ux_avg.append(tmp_ux[is_avg].mean())
        uy_avg.append(tmp_uy[is_avg].mean())
        ux_probe.append(tmp_ux[iProbe, jProbe])
        uy_probe.append(tmp_uy[iProbe, jProbe])
    tf_shore = dates.date2num(t_shore)

    # as Pandas data
    data_avg = pandas.DataFrame({'ux':ux_avg, 'uy':uy_avg}, index=t_shore)
    data_probe = pandas.DataFrame({'ux':ux_probe, 'uy':uy_probe}, index=t_shore)
    # resample
    data_avg = data_avg.resample('2S').mean().dropna()
    data_probe = data_probe.resample('2S').mean().dropna()

    # rotate
    if rotate:
        cosp = numpy.cos(numpy.deg2rad(rotate))
        sinp = numpy.sin(numpy.deg2rad(rotate))
        data_avg['vl'] = cosp*data_avg['ux'] - sinp*data_avg['uy']
        data_avg['vc'] = sinp*data_avg['ux'] + cosp*data_avg['uy']
        data_probe['vl'] = cosp*data_probe['ux'] - sinp*data_probe['uy']
        data_probe['vc'] = sinp*data_probe['ux'] + cosp*data_probe['uy']

    # spectra
    dt = 2. #must be the same as resample_short
    nwin = int(5*60*(1./dt)) #15 min x 60s/min x 1/dt sample/s
    offset = 10. #offset the cross-shore spectra
    f_vc_avg, Psd_vc_avg = signal.welch(data_avg['vc'], fs=1./dt, nperseg=nwin)
    f_vc_probe, Psd_vc_probe = signal.welch(data_probe['vc'], fs=1./dt, nperseg=nwin)

    #
    pyplot.figure()
    ax0 = pyplot.subplot(211, ylabel='cross (px)')
    pyplot.plot(data_avg.index, data_avg['vc'], 'k')
    pyplot.plot(data_probe.index, data_probe['vc'], 'r')
    pyplot.ylim(-1.5, 1.5)
    ax0.axhline(0., color='.2', ls='--', zorder=0)
    #
    ax1 = pyplot.subplot(212, ylabel='along (px)')
    pyplot.plot(data_avg.index, data_avg['vl'], 'r')
    pyplot.plot(data_probe.index, data_probe['vl'], 'k')
    pyplot.ylim(-1.5, 1.5)
    ax1.axhline(0., color='.2', ls='--', zorder=0)
    #
    fig = pyplot.figure()
    ax0 = fig.add_subplot(111, xlabel=r'f (Hz)', ylabel=r'PSD estimate',
                          xscale='log', yscale='log')
    ax0.plot(f_vc_avg[1:], Psd_vc_avg[1:], '-', color='k')
    ax0.plot(f_vc_probe[1:], Psd_vc_probe[1:], '-', color='r')

    pyplot.show()

if __name__=="__main__":
    #view('../tmp/test_v2/60m_jpg/20140613_0630_info.json')
    #test_cameradt()
    #test_misalignment()
    #compare_ADV_3()
    #ADVspectrum()
    #test_timeshift()
    #compare_estim()
    #compare_old()
    #export_series('../data')
    #test_median_vs_raw()
    #test_probe()
    #fulltest()
    #test_Dopplershift()
    test_Flow()
