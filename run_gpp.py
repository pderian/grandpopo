###
import glob
import os
import subprocess
###
from config import *
###

def process_group(pattern, params):
    """Extracts the frames from original videos and preprocess using decoder.py.

    Arguments:
        - pattern: the video search pattern, e.g. '20140312/12/*.mp4'.
        - params: dict of parameters for the preprocessing.

    Written by P. DERIAN 2017-04-05.
    """
    # check subdirectory according to filter options
    subdir = params['label']
    # and the parameters label
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
                              subdir, # according to parameters
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
                   f,
                   '-o', outdir,
                   '-p', prefix,
                   '-l', params['label'],
                   '-f', str(params['image_format']),
                   '-m', str(params['median_length']),
                   '-r', str(params['resolution']),
                   '-O', str(params['origin'][0]), str(params['origin'][1]),
                   '-d', str(params['dimensions'][0]), str(params['dimensions'][1]),
                   '-a', str(params['rotation']),
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
        # the area centered over sensors
        #pattern = '20140313/15/*.mp4' #sensors
        #process_group(pattern, params=PARAMS_COMP60)
        #process_group(pattern, params=PARAMS_COMP30)
        # the wide field for flash rip monitoring (TGRS paper)
        #pattern = '20140312/12/*.mp4' #flash rip
        #process_group(pattern, params=PARAMS_RIP120)
        # the wide field for swash (COASTAL DYN paper)
        pattern = '20140313/15/*.mp4' # [TMP]
        process_group(pattern, params=PARAMS_SWASH125)

    elif action=='estimator':
        pattern = '20140312/12/*.nc'
        #estim_group(pattern, which='raw')
        estim_group(pattern, which='median_5')
        #estim_group(pattern, which='radon_5')

