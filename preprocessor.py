###
import glob
import os
###
import numpy
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
import scipy.interpolate as interpolate
import scipy.io as scio
import scipy.signal as signal
import scipy.ndimage as ndimage
import scipy.spatial as spatial
import skimage.io as skio
import skimage.transform as sktransform
###
from config import *
###

### IMAGE PREPROCESSOR ###
##########################
class DataPreprocessor:
    """
    """
    def __init__(self, **kwargs):
        """
        Expected keywords:
            H: 3x3 homography (projection) matrix.
            origin: (x0, y0) the origin in [meter]
            dimensions: (xdim, ydim) the domain size in [meter]
            rotation: the rotation angle of the grid around (x0, y0) in [degree]
            resolution: the grid resolution in [meter].
            radon_length: the length of the radon filter, in [meter]. None or 0 skip the filtering

        Written by P. DERIAN 2016-01-23
        Modified by P. DERIAN 2017-01-11: added domain rotation
        """
        ### Store parameters
        self.param = kwargs
        # check radon parameters
        if ('radon_length' not in self.param) or (self.param['radon_length'] is None):
            self.param['radon_length'] = 0 # disable radon filtering
        self.param['radon_length_px'] = int(
            numpy.ceil(self.param['radon_length']/self.param['resolution']))
        # check median filter parameters
        if ('median_length' not in self.param) or (self.param['median_length'] is None):
            self.param['median_length'] = 0 # disable radon filtering
        self.param['median_length_px'] = int(
            numpy.ceil(self.param['median_length']/self.param['resolution']))
        ### Load the projection matrix and create the corresponding transform
        self.H = numpy.array(self.param['H'])
        self.projection = sktransform.ProjectiveTransform(self.H)
        ### Real world interpolation coordinates
        self.X, self.Y = domain_grid(self.param['origin'], self.param['dimensions'],
                                     self.param['rotation'], self.param['resolution'])
        YX = numpy.hstack((self.Y.ravel().reshape((-1,1)),
                           self.X.ravel().reshape((-1,1)),
                           )) # interpolation coordinates
        self.shape = self.X.shape # this is the grid shape
        ### Image interpolation coordinates
        self.iYX = self.projection.inverse(YX)

    def __call__(self, args):
        """
        For use with multiprocessing.

        Written by P. DERIAN 2016-05-20
        """
        return self.process_image(*args)

    def process_image(self, image, as_uint=False):
        """
        Main processing pipe.

        Arguments:
            - image: a (M,N) (graysacle) or (M,N,3) (RGB) image;
            - as_uint=False: convert output to valid [0, 255] uint8 images.
        Return: (img_rect, img_filt)
            - img_rect: the rectified image
            - img_filt: the filtered version of img_rect.

        Written by P. DERIAN 2016-03-09
        Modified by P. DERIAN 2017-01-05: added clipping and int conversion.
        """
        # first rectify and grid
        img_rect = self.grid_image(image)
        img_filt = None # the default filtered version
        # then apply Radon
        if self.param['radon_length_px']>1:
            img_filt = radon_filter(
                img_rect,
                self.param['radon_length_px'],
                angles=numpy.linspace(0., 180., num=img_processed.shape[0], endpoint=False),
                )
        # or median fiter
        elif self.param['median_length_px']>1:
            img_filt = ndimage.median_filter(
                img_rect,
                size=self.param['median_length_px'],
                mode='reflect',
                ) # this is the low-pass filter
        # if as_uint: clip and convert back to uint8
        if as_uint:
            numpy.clip(img_rect, 0., 1., img_rect)
            img_rect = (255.*img_rect).astype('uint8')
            if img_filt is not None:
                numpy.clip(img_filt, 0., 1., img_filt)
                img_filt = (255.*img_filt).astype('uint8')
        return img_rect, img_filt

    def grid_image(self, image):
        """
        Grid the supplied image on the domain.

        Arguments:
            image a 2D grayscale or RBG image.
        Output:
            the 2D re-gridded image.

        Written by P. DERIAN 2016-01-23
        """
        # interpolate, reshape and return
        if image.ndim==2:
            return ndimage.interpolation.map_coordinates(image, self.iYX.T).reshape(self.shape)
        # TODO: make it faster???
        elif image.ndim==3:
            return numpy.dstack(
                (ndimage.interpolation.map_coordinates(image[:,:,i], self.iYX.T).reshape(self.shape)
                 for i in xrange(3)),
                )

    def demo(self, imagefile, outfile=None):
        """
        Show the preprocessor output.

        Arguments:
            imagefile: path to an image to be processed
        """
        def get_yx_polygon(X, Y):
            return numpy.array([[Y[iy, ix], X[iy, ix]]
                               for [iy, ix] in zip([0,-1,-1,0], [0,0,-1,-1])])

        # load image
        img = skio.imread(imagefile, as_grey=True)
        # interpolate
        imgc = self.grid_image(img)

        # create figure
        dpi = 90.
        fig, (ax1, ax2) = pyplot.subplots(2,1, figsize=(1000./dpi, 1000./dpi))

        ### interpolation boundaries
        yxBound = get_yx_polygon(self.X, self.Y)
        iyxBound = self.projection.inverse(yxBound)

        ### compas
        x0Compas = self.X[0,0] + 110. #[m]
        y0Compas = self.Y[0,0] + 80. #[m]
        lCompas = 10. #[m]
        x1Compas = x0Compas
        y1Compas = y0Compas + lCompas
        x2Compas = x0Compas + lCompas
        y2Compas = y0Compas
        yxCompas = numpy.array([[y0Compas, x0Compas],
                    [y1Compas, x1Compas],
                    [y2Compas, x2Compas]])
        iyxCompas = self.projection.inverse(yxCompas)
        iyxCompas -= iyxCompas[0] - [[550, 1000]] #shift the arrows
        yxCompas -= yxCompas[0] - [[self.Y[0,0]+5., self.X[0,0]+5]]
        for i in [1,2]:
            iyxCompas[i,:] -= iyxCompas[0,:]
            iyxCompas[i,:] *= 30./numpy.sqrt(numpy.sum(iyxCompas[i,:]**2))
        print iyxCompas

        ### aquapro & ADV sensors
        iyxAquapro = numpy.array([[370, 600],]) #pixels
        yxAquapro = self.projection(iyxAquapro)
        iyxADV = numpy.array([[360, 660],]) #pixels
        yxADV = self.projection(iyxADV)
        iyxRelease = numpy.array([[270, 1175],]) #pixels
        yxRelease = self.projection(iyxRelease)

        ### subdomain
        # 60-m estimation area around sensors
        X60, Y60 = domain_grid(PARAMS_COMP60['origin'], PARAMS_COMP60['dimensions'],
                               PARAMS_COMP60['rotation'], PARAMS_COMP60['resolution'])
        yx60 = get_yx_polygon(X60, Y60)
        iyx60 = self.projection.inverse(yx60) # pixel coords
        y60_label = yx60[:,0].min() # world coord for the text label
        x60_label = yx60[:,1].mean()
        area60_color = 'purple'
        # 30-m area around sensors
        X30, Y30 = domain_grid(PARAMS_COMP30['origin'], PARAMS_COMP30['dimensions'],
                               PARAMS_COMP30['rotation'], PARAMS_COMP30['resolution'])
        yx30 = get_yx_polygon(X30, Y30)
        iyx30 = self.projection.inverse(yx30) # pixel coords
        y30_label = yx30[:,0].min() # world coord for the text label
        x30_label = yx30[:,1].mean()
        area30_color = 'orange'
        # 120x60 m wide area for flash rip
        X120, Y120 = domain_grid(PARAMS_RIP120['origin'], PARAMS_RIP120['dimensions'],
                                 PARAMS_RIP120['rotation'], PARAMS_RIP120['resolution'])
        yx120 = get_yx_polygon(X120, Y120)
        iyx120 = self.projection.inverse(yx120) # pixel coords
        y120_label = yx120[:,0].min() # world coord for the text label
        x120_label = yx120[:,1].mean()
        area120_color = 'c'
        # 125x45 m for COASTDYN
        X125, Y125 = domain_grid(PARAMS_SWASH125['origin'], PARAMS_SWASH125['dimensions'],
                                 PARAMS_SWASH125['rotation'], PARAMS_SWASH125['resolution'])
        yx125 = get_yx_polygon(X125, Y125)
        iyx125 = self.projection.inverse(yx125) # pixel coords
        y125_label = Y125[0,:].mean() # world coord for the text label
        x125_label = yx125[:,1].mean()
        area125_color = 'pink'



        #### original data
        ax1.set_title('Original')
        ax1.set_xlabel('i [px]')
        ax1.set_ylabel('j [px]')
        ax1.imshow(img, cmap='gray', interpolation='nearest', vmin=0.1, vmax=0.9)
        # plot the various areas and add a marker at the origin
        # Note: roll to get x,y instead of y, x
        ax1.add_artist(patches.Polygon(numpy.roll(iyxBound,1,axis=-1), fill=False, ls='--'))
        ax1.add_artist(patches.Polygon(numpy.roll(iyx60,1,axis=-1), fill=False, color=area60_color))
        ax1.plot(iyx60[0,1], iyx60[0,0], 'o', color=area60_color)
        ax1.add_artist(patches.Polygon(numpy.roll(iyx30,1,axis=-1), fill=False, color=area30_color))
        ax1.plot(iyx30[0,1], iyx30[0,0], 'o', color=area30_color)
        ax1.add_artist(patches.Polygon(numpy.roll(iyx120,1,axis=-1), fill=False, color=area120_color))
        ax1.plot(iyx120[0,1], iyx120[0,0], 'o', color=area120_color)
        ax1.add_artist(patches.Polygon(numpy.roll(iyx125,1,axis=-1), fill=False, color=area125_color))
        ax1.plot(iyx125[0,1], iyx125[0,0], 'o', color=area125_color)
        # the sensors
        ax1.plot(iyxAquapro[0,1], iyxAquapro[0,0], '*r')
        ax1.plot(iyxADV[0,1], iyxADV[0,0], '*g')
        ax1.plot(iyxRelease[0,1], iyxRelease[0,0], '*b')
        # the compas
        ax1.quiver([iyxCompas[0,1],],
                   [iyxCompas[0,0],],
                   [iyxCompas[1,1],],
                   [iyxCompas[1,0],],
                   color=['w', 'k'],
                   units='xy', angles='xy',
                   scale_units='xy',
                   )
        ax1.text(iyxCompas[0,1]-5, iyxCompas[0,0]-5, 'N', color='w', weight='bold',
                 va='bottom', ha='right')
        #
        ax1.set_xlim(0, img.shape[1])
        ax1.set_ylim(img.shape[0], 0)
        ### rectified data
        ax2.set_aspect('equal')
        ax2.set_title('Rectified, {} m/px'.format(self.param['resolution']))
        ax2.set_xlabel('easting x [m]')
        ax2.set_ylabel('northing y [m]')
        ax2.pcolormesh(self.X, self.Y, imgc, cmap='gray')
        # the areas
        ax2.add_artist(patches.Polygon(numpy.roll(yx60,1,axis=-1), fill=False, color=area60_color))
        ax2.plot(yx60[0,1], yx60[0,0], 'o', color=area60_color)
        ax2.text(x60_label, y60_label, 'TGRS(v1) {}x{} m'.format(*PARAMS_COMP60['dimensions']),
                va='top', ha='center', color=area60_color)
        ax2.add_artist(patches.Polygon(numpy.roll(yx30,1,axis=-1), fill=False, color=area30_color))
        ax2.plot(yx30[0,1], yx30[0,0], 'o', color=area30_color)
        ax2.text(x30_label, y30_label, 'TGRS(v2?) {}x{} m'.format(*PARAMS_COMP30['dimensions']),
                 va='top', ha='center', color=area30_color)
        ax2.add_artist(patches.Polygon(numpy.roll(yx120,1,axis=-1), fill=False, color=area120_color))
        ax2.plot(yx120[0,1], yx120[0,0], 'o', color=area120_color)
        ax2.text(x120_label, y120_label, 'TGRS(v1,v2) {}x{} m'.format(*PARAMS_RIP120['dimensions']),
                 va='top', ha='center', color=area120_color)
        ax2.add_artist(patches.Polygon(numpy.roll(yx125,1,axis=-1), fill=False, color=area125_color))
        ax2.plot(yx125[0,1], yx125[0,0], 'o', color=area125_color)
        ax2.text(x125_label, y125_label, 'COASTDYN {}x{} m'.format(*PARAMS_SWASH125['dimensions']),
                 va='top', ha='center', color=area125_color)
        ax2.add_artist(
            patches.Circle((AVG_PROBE['x'], AVG_PROBE['y']), radius=AVG_PROBE['r'],
                           fill=True, color='teal', alpha=0.75),
            )
        # the sensors
        ax2.plot(yxAquapro[0,1], yxAquapro[0,0], '*r')
        ax2.text(yxAquapro[0,1]+1, yxAquapro[0,0]+1, 'Aquapro', va='baseline', ha='left', color='r')
        ax2.plot(yxADV[0,1], yxADV[0,0], '*g')
        ax2.text(yxADV[0,1]-1, yxADV[0,0]-1, 'ADV', va='baseline', ha='right', color='g')
        ax2.plot(yxRelease[0,1], yxRelease[0,0], '*b')
        ax2.text(yxRelease[0,1], yxRelease[0,0]-1.5, 'Dye release', va='baseline', ha='center', color='b')
        # the compas
        ax2.quiver([yxCompas[0,1],],
                   [yxCompas[0,0],],
                   [yxCompas[1,1]-yxCompas[0,1],],
                   [yxCompas[1,0]-yxCompas[0,0],],
                   color=['w', 'k'],
                   units='xy', angles='xy', scale_units='xy')
        ax2.text(yxCompas[0,1], yxCompas[0,0]-1., 'N', color='w', weight='bold',
                 va='top', ha='center')
        #
        ax2.set_xlim(yxBound[:,1].min(), yxBound[:,1].max())
        ax2.set_ylim(yxBound[:,0].min(), yxBound[:,0].max())

        print 'Aquapro (x,y) = ({:.0f}, {:.0f}) (px) -> ({:.2f}, {:.2f}) (m)'.format(
            iyxAquapro[0,1], iyxAquapro[0,0], yxAquapro[0,1], yxAquapro[0,0])
        print 'ADV (x,y) = ({:.0f}, {:.0f}) (px) -> ({:.2f}, {:.2f}) (m)'.format(
            iyxADV[0,1], iyxADV[0,0], yxADV[0,1], yxADV[0,0])

        pyplot.subplots_adjust(left=.07, bottom=.05, right=.97, top=.96)
        if outfile is not None:
            fig.savefig(outfile, dpi=dpi)
            print 'saved', outfile
        pyplot.show()

### HELPERS ###
###############

def radon_filter(image, length, angles=numpy.arange(180)):
    """
    High-pass filter an image to remove structures bigger than "length" pixels,
    using the radon transform.

    Arguments:
        image: a 2D image.
        length: the filter length in [pixel].
        angles: array of angles in [degree] for the radon transform.
                Default to numpy.arange(180).
    Output:
        the filtered image, SQUARE!.

    Written 2016-01-22 by P. DERIAN - www.pierrederian.net
    """
    # Remove the image mean (basic detrend)
    image = image - image.mean()
    # Apply radon transform
    # Note: rotation axis is at sinogram.shape[0]//2 along the first dimension (see doc)
    # Note: the first dimension is numpy.ceil(numpy.sqrt(2.)*max(image.shape)) (i.e. longest diagonal)
    sinogram = sktransform.radon(image, theta=angles)
    ############################
    ### compute the envelope ###
    # Note: USELESS in the end. Left here for posterity since it wasn't so trivial.
    if 0:
        # The transform use a padded square image, so we use the longest original image size.
        m = max(image.shape)
        # The origin of the transform lies at the center.
        xc = m//2 # in [pixel]
        # For any given projection angle, the envelop coordinate is given by the max distance between
        # - a point of the image,
        # - and its projection on the line that passes through the origin,
        #   forming said angle (increasing counterclockwise) w.r.t. the horizontal.
        # The furthest points are always the image corners,
        # but because of the central symmetry that applies to this situation,
        # we only check the upper half corners here.
        corners = [(m/2.,m/2.), (-m/2.,m/2.)]
        # Slope s of the line y=sx at every projection angle.
        angles_rad = numpy.deg2rad(angles) # in [radian]
        angles_sin = numpy.sin(angles_rad)
        angles_cos = numpy.cos(angles_rad)
        slope = numpy.where(numpy.abs(angles_sin)<numpy.abs(angles_cos), angles_sin, angles_cos)
        # For every slope (i.e. every projection line), we compute the distance
        # from the corner to the corresponding projection line
        # and take the max.
        # https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
        envelop = numpy.array(
            [max([numpy.abs(s*c[0]-c[1])/(numpy.sqrt(1.+s**2)) for c in corners]) for s in slope],
            )
    ############################
    # now high-pass filter the sinogram along columns
    sinogram -= ndimage.uniform_filter1d(sinogram, length, axis=0, mode='constant')
    # and rebuild the image
    image_filtered = sktransform.iradon(sinogram, output_size=max(image.shape),
                                        interpolation='linear', filter='hann')
    return image_filtered

def estimate_H_from_mapping(matfile, sstep=100):
    """
    Estimate the projective transform from the coordinates in the .mat file.
    Uses the pixel->world mapping defined in the mat file
    to recover the projective transform.

    Arguments:
        matfile: the "CoordonneesXY_GPP_20141013.mat" Matlab file
                 containing 'X' and 'Y' arrays of world coordinates.
        sstep: sub-sampling step applied to points
               before estimating the transform (for speed-up).
    Output:
        the 3x3 projection matrix

    Written by P. DERIAN 2016-01-23
    """
    pmat = scio.loadmat(matfile)
    X = pmat['X']
    Y = pmat['Y']
    # pixel coords
    ix = numpy.arange(pmat['X'].shape[1])
    iy = numpy.arange(pmat['X'].shape[0])
    iX, iY = numpy.meshgrid(ix, iy)
    # downsamling step
    iYX = numpy.hstack((iY[::sstep, ::sstep].ravel().reshape((-1,1)),
                        iX[::sstep, ::sstep].ravel().reshape((-1,1))))
    YX = numpy.hstack((Y[::sstep, ::sstep].ravel().reshape((-1,1)),
                       X[::sstep, ::sstep].ravel().reshape((-1,1))))
    # estimate the transform
    projection = sktransform.ProjectiveTransform()
    if projection.estimate(iYX, YX):
        return projection.params
    else:
        print '[!] projection estimation failed'
        return None

def check_H(matfile, H, sstep=100):
    """
    Check the accuracy of the projective transform H.
    Uses the pixel->world mapping defined in the mat file
    to verify the transform.

    Arguments:
        H: the projection matrix.
        matfile: the "CoordonneesXY_GPP_20141013.mat" Matlab file
                 containing 'X' and 'Y' arrays of world coordinates.
        sstep: sub-sampling step applied to points (for speed-up).

    Return: errProj, errInv the projection and inverse projection RMS errors.

    Note: RMS errors on the order of 1e-10 are found with H estimated using a sstep=100.
        This ought to be be accurate enough...

    Written by P. DERIAN 2016-03-23
    """

    # create the projection
    projection = sktransform.ProjectiveTransform(H)

    # input points
    pmat = scio.loadmat(matfile)
    ix = numpy.arange(pmat['X'].shape[1]) #pixels
    iy = numpy.arange(pmat['X'].shape[0])
    iX, iY = numpy.meshgrid(ix, iy)
    tX = pmat['X'] # true world (UTM coords)
    tY = pmat['Y']
    # downsamling step
    iYX = numpy.hstack((iY[::sstep, ::sstep].ravel().reshape((-1,1)),
                        iX[::sstep, ::sstep].ravel().reshape((-1,1))))

    tYX = numpy.hstack((tY[::sstep, ::sstep].ravel().reshape((-1,1)),
                        tX[::sstep, ::sstep].ravel().reshape((-1,1))))
    # compute projection of pixels
    p_iYX = projection(iYX)
    errProj = numpy.sqrt(numpy.mean(numpy.sum(numpy.square(tYX - p_iYX), axis=-1))) # RMS error (meters)
    # compute inverse projection of UTM
    ip_tYX = projection.inverse(tYX)
    errInv = numpy.sqrt(numpy.mean(numpy.sum(numpy.square(iYX - ip_tYX), axis=-1))) # RMS error (pixels)

    print 'projection RMS error (meters):', errProj
    print 'inverse proj. RMS error (pixels):', errInv
    return errProj, errInv

def test_projection(matfile, H, imgfile):
    """ another test to debug projection.
    This one alos uses an image for comparison

    Written by P. DERIAN 2016-05-31
    """
    # mat data
    pmat = scio.loadmat(matfile)
    ix = numpy.arange(pmat['X'].shape[1]) #pixels
    iy = numpy.arange(pmat['X'].shape[0])
    tx = pmat['X'] #UTM
    ty = pmat['Y']
    # homemade preprocessor
    preprocessor = DataPreprocessor(
        H=DEFAULT_H,
        xbounds=(370120., 370380.),
        ybounds=(694000., 694150.),
        resolution=0.2,
        )
    # load image
    img = skio.imread(imgfile)
    # grid with preprocessor
    imc = preprocessor.grid_image(img)
    skio.imsave("../tmp/preproc_img.jpg", imc)
    # now from the mat data (nearest neighbor)
    tXY = numpy.hstack((ty.ravel().reshape((-1,1)),
                        tx.ravel().reshape((-1,1))))
    X,Y = numpy.meshgrid(preprocessor.x, preprocessor.y)
    XY = numpy.hstack((Y.ravel().reshape((-1,1)),
                       X.ravel().reshape((-1,1))))
    tree = spatial.cKDTree(tXY)
    _, iinterp = tree.query(XY)
    for i in xrange(3):
        tmp = img[:,:,i].flat
        imc[:,:,i] = tmp[iinterp].reshape((imc.shape[0], imc.shape[1]))
    skio.imsave("../tmp/nninterp_img.jpg", imc)

def test_mapping(matfile, H):
    """
    Test the mapping provided in the reference matfile.

    Written by P. DERIAN 2016-03-23
    """

    pmat = scio.loadmat(matfile)
    ix = numpy.arange(pmat['X'].shape[1]) #pixels
    iy = numpy.arange(pmat['X'].shape[0])
    iX, iY = numpy.meshgrid(ix, iy)
    tX = pmat['X'] # true world (UTM coords)
    tY = pmat['Y']

    ref_ix, ref_iy = 248, 398
    ref_x, ref_y = 370359.89, 694046.30

    projection = sktransform.ProjectiveTransform(H)
    p_yx = projection([[ref_iy, ref_ix],])
    print p_yx[0,1], p_yx[0,0]

    print 'XHorPixAqua: {} -> to UTM: {:.2f} (XUTMAqua: {:.2f})'.format(ref_ix, tX[ref_iy, ref_ix], ref_x)
    print 'YVertPixAqua: {} -> to UTM: {:.2f} (YUTMAqua: {:.2f})'.format(ref_iy, tY[ref_iy, ref_ix], ref_y)

### MAIN ####
#############
if __name__=="__main__":

    ### tests
    #H = DEFAULT_H
    #mappingFile = '../data/GrandPopo_Misc/CoordonneesXY_GPP_20141013.mat'
    #H = estimate_H_from_mapping(mappingFile)
    #print numpy.array_repr(H, precision=16)
    #check_H(mappingFile, H)
    #test_mapping(mappingFile, H)
    #test_projection(mappingFile, H, '../data/GrandPopo_Misc/ex_rgb_frame.jpg')

    preprocessor = DataPreprocessor(
        H=DEFAULT_H,
        origin=(370220., 694040.),
        dimensions=(150., 75.),
        rotation=0.,
        resolution=0.2,
        )
    preprocessor.demo("resources/sample_frame_release.jpg",
                      "resources/grandpopo_config.jpg")
