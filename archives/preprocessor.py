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

### this H was estimated from the mapping in the mat file.
# see estimate_H_from_mapping()
DEFAULT_H = numpy.array(
    [[  2.4795542480487102e+03,   6.7539967752787433e+01,
       -2.2440068592504461e+05],
     [  1.3228681378707056e+03,   3.5955050912544941e+01,
       -1.1965298848769323e+05],
     [  3.5719691693753153e-03,   9.7256605585734505e-05,
       -3.2313166347066136e-01]],
    )

### IMAGE PREPROCESSOR ###
##########################
class DataPreprocessor:
    """
    """
    def __init__(self, **kwargs):
        """
        Expected keywords:
            H: 3x3 homography (projection) matrix.
            xbounds: (xmin, xmax) the x (horizontal) domain of interest in [meter].
            ybounds: (ymin, ymax) the y (vertical) domain of interest in [meter].
            resolution: the grid resolution in [meter].
            radon_length: the length of the radon filter, in [meter]. None or 0 skip the filtering

        Written by P. DERIAN 2016-01-23
        Modified by P. DERIAN 2016-03-09: added radon transform capability
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
        self.H = self.param['H']
        self.projection = sktransform.ProjectiveTransform(self.H)
        ### Real world interpolation coordinates
        self.x = numpy.arange(self.param['xbounds'][0], self.param['xbounds'][1],
                               step=self.param['resolution'],
                               ) # 1d
        self.y = numpy.arange(self.param['ybounds'][0], self.param['ybounds'][1],
                               step=self.param['resolution'],
                               ) # 1d
        self.shape = (self.y.size, self.x.size) # shape of the gridded data
        X, Y = numpy.meshgrid(self.x, self.y) # 2d
        YX = numpy.hstack((Y.ravel().reshape((-1,1)),
                           X.ravel().reshape((-1,1)),
                           )) # interpolation coordinates
        ### Image interpolation coordinates
        self.iYX = self.projection.inverse(YX)

    def __call__(self, args):
        """
        For use with multiprocessing.

        Written by P. DERIAN 2016-05-20
        """
        return self.process_image(*args)

    def process_image(self, image):
        """
        Main processing pipe.

        Written by P. DERIAN 2016-03-09
        """
        # first grid
        img_processed = self.grid_image(image)
        # then apply Radon
        if self.param['radon_length_px']>1:
            print 'radon'
            img_processed = radon_filter(
                img_processed,
                self.param['radon_length_px'],
                angles=numpy.linspace(0., 180., num=img_processed.shape[0], endpoint=False),
                )
        # or median fiter
        elif self.param['median_length_px']>1:
            img_processed -= img_processed.mean()
            img_med = ndimage.median_filter(
                img_processed,
                size=self.param['median_length_px'],
                mode='reflect',
                ) # this is the low-pass filter
            img_processed -= img_med # now we retrieve the high-pass filter
        # ...
        return img_processed

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

    def demo(self, imagefile):
        """
        Show the preprocessor output.

        Arguments:
            imagefile: path to an image to be processed
        """
        # load image
        img = skio.imread(imagefile, as_grey=True)
        imgc = self.grid_image(img)
        # interpolate
        fig, (ax1, ax2) = pyplot.subplots(1,2)
        # interpolation boundaries
        xmin, xmax = self.param['xbounds']
        ymin, ymax = self.param['ybounds']
        yxBound = numpy.array([[ymax, xmin], [ymax, xmax], [ymin, xmax], [ymin, xmin]])
        iyxBound = self.projection.inverse(yxBound)
        # aquapro
        iyxAquapro = numpy.array([[370, 600],]) #pixels
        yxAquapro = self.projection(iyxAquapro)
        iyxADV = numpy.array([[360, 660],]) #pixels
        yxADV = self.projection(iyxADV)
        # display
        ax1.set_title('Original')
        ax1.imshow(img, cmap='gray', interpolation='none')
        # roll to get x,y instead of y, x
        ax1.add_artist(patches.Polygon(numpy.roll(iyxBound,1,axis=-1), fill=False))
        ax1.plot(iyxAquapro[0,1], iyxAquapro[0,0], '*r')
        ax1.plot(iyxADV[0,1], iyxADV[0,0], '*g')
        ax1.set_xlim(0, img.shape[1])
        ax1.set_ylim(img.shape[0], 0)
        #
        ax2.set_title('Interpolated, {} m/px'.format(self.param['resolution']))
        ax2.imshow(imgc, cmap='gray', interpolation='none',
            extent=[self.x[0], self.x[-1], self.y[-1], self.y[0]])
        ax2.plot(yxAquapro[0,1], yxAquapro[0,0], '*r')
        ax2.plot(yxADV[0,1], yxADV[0,0], '*g')
        ax2.set_xlim(self.x[0], self.x[-1])
        ax2.set_ylim(self.y[-1], self.y[0])

        print 'Aquapro (x,y) = ({:.0f}, {:.0f}) (px) -> ({:.2f}, {:.2f}) (m)'.format(
            iyxAquapro[0,1], iyxAquapro[0,0], yxAquapro[0,1], yxAquapro[0,0])
        print 'ADV (x,y) = ({:.0f}, {:.0f}) (px) -> ({:.2f}, {:.2f}) (m)'.format(
            iyxADV[0,1], iyxADV[0,0], yxADV[0,1], yxADV[0,0])

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
    H = DEFAULT_H
    mappingFile = '../data/GrandPopo_Misc/CoordonneesXY_GPP_20141013.mat'
    #H = estimate_H_from_mapping(mappingFile)
    #print numpy.array_repr(H, precision=16)
    #check_H(mappingFile, H)
    #test_mapping(mappingFile, H)
    test_projection(mappingFile, H, '../data/GrandPopo_Misc/ex_rgb_frame.jpg')

    preprocessor = DataPreprocessor(
        H=DEFAULT_H,
        xbounds=(370250., 370350.),
        ybounds=(694050., 694130.),
        resolution=0.1,
        )
    #preprocessor.demo("tmp/raw_000.png")
