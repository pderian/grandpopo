###
import glob
import os
import subprocess
###
ROOT_RAWDATA_DIR = '/Volumes/LaCie_Mac/pderian/data_GPP/Videos'
ROOT_PREPROC_DIR = '/Volumes/LaCie_Mac/pderian/data_GPP/Frames'
###
# Parameters
#x_UTM = 370325 #sensors
#y_UTM = 694098
#x_UTM = 370260 #flash rip
#y_UTM = 694090
#dim = 60.
#resolution = .1 #sensors
#resolution = .2 #flash rip
#xmin, xmax = x_UTM - dim/2., x_UTM + dim/2. #sensors
#xmin, xmax = x_UTM - 20., x_UTM + 100. #flash rip
#ymin, ymax = y_UTM - 3*dim/4., y_UTM + dim/4.


def process_group(pattern, median=0.):
    """Extracts the frames from original videos and preprocess using decoder.py.

    Arguments:
        - pattern: the video search pattern, e.g. '20140312/12/*.mp4'.
        - median=0.: the length in [m] of the median filter to apply.

    [TODO] add bounds & resolution options.

    Written by P. DERIAN 2017-04-05.
    """
    # check subdirectory according to filter options
    subdir = 'median_{}m'.format(int(median)) if median>0 else 'raw'
    print 'Processing:', subdir

    # search for videos matching the pattern
    search = os.path.join(ROOT_RAWDATA_DIR, pattern)
    print 'Search pattern:', search
    flist = sorted(glob.glob(search))

    # for each matching video
    for f in flist:
        # video structures (copied from LEGOS FTP) is yyyymmdd/HH/MM.mp4
        # and we want to store frames as yyyymmdd/yyyymmdd_HH/yyyymmdd_HHMM/yyyymmdd_HHMM_<index>.<format>
        # so: recursively split to extract basename, hour and date
        p, fname = os.path.split(f)
        p, hour = os.path.split(p)
        p, date = os.path.split(p)
        minute, _ = os.path.splitext(fname)
         # compute output dir, and prefix for frames
        outdir = os.path.join(ROOT_PREPROC_DIR,
                              subdir, # according to filter option
                              date,
                              '{}_{}'.format(date, hour),
                              '{}_{}{}'.format(date, hour, minute),
                              )
        prefix = '{}_{}{}_'.format(date, hour, minute)
        # create output directory if neeeded
        if not os.path.exists(outdir):
            print 'Creating output directory', outdir
            os.makedirs(outdir, 0755)
        # call decoder
        command = ['python', '-u', 'decoder.py',
                   '--img',
                   '-o', outdir,
                   '-p', prefix,
                   '-f', 'jpg',
                   '-m', str(median),
                   f,
                   ]
        subprocess.call(command)

def estim_group(pattern, which='raw', serverPort=9999):

    root_preproc_dir = '../data/GrandPopo_Preprocessed'
    root_estim_dir = '../estim/GrandPopo'

    # check parameter
    if which not in ['raw', 'radon_5', 'median_5']:
        raise RuntimeError("unknown dataset: {}".format(which))

    # find data files
    flist = sorted(glob.glob(os.path.join(root_preproc_dir, which, pattern)))
    print flist

    # get their relative path
    outdirs = []
    for f in flist:
        # recursively split to extract basename, hour and date
        p, fname = os.path.split(f)
        p, hour = os.path.split(p)
        p, date = os.path.split(p)
        # compute output dir, nc file
        outdirs.append(os.path.join(root_estim_dir, which, date, hour))
    # remove duplicates
    subsets = set(outdirs)
    # tupes of (input files, output dirs)
    io = zip(flist, outdirs)

    # for each subset (= unique output directory)
    for sub in subsets:
        # create output directory as needed
        if not os.path.exists(sub):
            print 'Creating output directory', sub
            os.makedirs(sub, 0700)
        # select files for this subset
        subfiles = [i for i,o in io if o==sub]
        # call estimator
        # typically: python estimator.py 9999 -f ../data/GrandPopo_Preprocessed/raw/20140313/06/30.nc -o ../tmp -c 'test estim raw'
        command = ['python', 'estimator.py',
                   str(serverPort),
                   '-o', sub,
                   '-c', 'GrandPopo surface current estimation - "{}" data.'.format(which),
                   '-f',] + subfiles
        subprocess.call(command)


if __name__=="__main__":

    #action = 'estimator'
    action = 'decoder'

    if action=='decoder':
        #pattern = '20140312/12/*.mp4' #flash rip
        pattern = '20140313/15/00.mp4' #sensors
        process_group(pattern, median=5)
    elif action=='estimator':
        pattern = '20140312/12/*.nc'
        #estim_group(pattern, which='raw')
        estim_group(pattern, which='median_5')
        #estim_group(pattern, which='radon_5')

