###
import argparse
import datetime
import glob
import json
import os
import sys
###
import matplotlib.pyplot as pyplot
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



if __name__=="__main__":
    view('../tmp/test_v2/60m_jpg/20140613_0630_info.json')
