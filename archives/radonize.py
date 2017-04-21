import argparse
import sys
###
import numpy
import skimage.io as skio
import skimage.exposure as skexp
###
from preprocessor import radon_filter
###

def main(argv):
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", help="image input file")
    parser.add_argument("-o", "--output", default=None, dest="outputfile", help="image output file")
    parser.add_argument("-l", "--length", default=10, type=int, dest="radon_length", help="length of radon filter in [px]")
    parser.add_argument("-c", "--channel", dest="channel", default=None, type=int, choices=[0, 1, 2], help="channel to be used")
    parser.add_argument("-n", "--nangles", dest="nangles", default=None, type=int, help="number of angles")
    parser.add_argument("--crop", dest="crop", nargs=2, default=None, type=int, help="'width height' cropping dimensions")
    parser.add_argument("--offset", dest="offset", nargs=2, default=[0, 0], type=int, help="'offset_width offset_height' cropping offset")
    parser.add_argument("--scale", dest="scale", nargs=2, default=None, type=float, help="'scale_min scale_max' scaling parameters")
    parser.add_argument("-s", "--show", dest="show", action="store_true", help="show images")
    parser.set_defaults(show=False)
    args = parser.parse_args(argv)

    # load image
    img_raw = skio.imread(args.inputfile)
    # maybe crop
    if args.crop:
        img_raw = img_raw[args.offset[1]:args.offset[1]+args.crop[1], args.offset[0]:args.offset[0]+args.crop[0], :]
    # convert to grayscale    
    if args.channel is None:
        # usual luminance formula
        weights = [0.2126, 0.7152, .0722]
        if img_raw.ndim==4:
            weights.append(1.);
        img_in = numpy.average(img_raw, axis=-1, weights=weights)
    else:
        # use specified channel
        img_in = img_raw[:,:,args.channel]
    # rescale to 0, 1
    img_in = skexp.rescale_intensity(img_in.astype(float), in_range=(0., 255.), out_range=(0., 1.))
    # apply radon filter
    if args.nangles is None:
        args.nangles = min(img_in.shape)
    img_out = radon_filter(img_in, args.radon_length,
                           angles=numpy.linspace(0., 180, args.nangles, endpoint=False))
    # rescale
    if args.scale:
        img_out = skexp.rescale_intensity(img_out, in_range=tuple(args.scale), out_range=(0., 255.))
    # show
    if args.show:
        import matplotlib.pyplot as pyplot
        pyplot.subplot(131)
        pyplot.imshow(img_raw)
        pyplot.subplot(132)
        p = pyplot.imshow(img_in, cmap='gray')
        pyplot.colorbar(p)
        pyplot.subplot(133)
        p = pyplot.imshow(img_out, cmap='gray')
        pyplot.colorbar(p)
        pyplot.show()
    # export
    if args.outputfile:
        skio.imsave(args.outputfile, img_out.astype(int))
        print 'saved', args.outputfile
    

if __name__=="__main__":
    main(sys.argv[1:])
