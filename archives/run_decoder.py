###
import glob
import os
import subprocess
###

def process_group(pattern, radon=0., median=0., force_replace=False):
    root_rawdata_dir = '../data/GrandPopo_Videos'
    root_preproc_dir = '../data/GrandPopo_Preprocessed'

    # check parameters
    if not (radon or median):
        # raw (no filtering)
        subdir = 'raw'
    elif radon and (not median):
        # radon filtering
        subdir = 'radon_{}'.format(int(radon))
    elif median and (not radon):
        # median filtering
        subdir = 'median_{}'.format(int(median))
    else:
        # error
        raise RuntimeError("radon and median filtering cannot be ebabled simultaneously")
    print 'Processing:', subdir 

    search = os.path.join(root_rawdata_dir, pattern)
    print search
    flist = sorted(glob.glob(search))
    for f in flist:
        # recursively split to extract basename, hour and date
        p, fname = os.path.split(f)
        p, hour = os.path.split(p)
        p, date = os.path.split(p)
        bname, _ = os.path.splitext(fname)
         # compute output dir, nc file
        outdir = os.path.join(root_preproc_dir, subdir, date, hour)
        outfile = os.path.join(outdir, bname+'.nc')
        # create output directory if neeeded
        if not os.path.exists(outdir):
            print 'Creating output directory', outdir
            os.makedirs(outdir, 0700)
        # check if file already exists
        if os.path.isfile(outfile):
            if force_replace:
                print 'overwriting', outfile
            else:
                #raise RuntimeError('file {} already exists, aborting.'.format(outfile))
                print 'file {} exists, skipping.'.format(outfile)
                continue
        # call decoder
        command = ['python', '-u', 'decoder.py',
                   '-o', outfile,
                   '-m', str(median),
                   '-r', str(radon),
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
        replace = True
        #process_group(pattern, radon=0, median=0, force_replace=replace)
        process_group(pattern, radon=5, median=0, force_replace=replace)
        #process_group(pattern, radon=0, median=5, force_replace=replace)
    elif action=='estimator':
        pattern = '20140312/12/*.nc'
        #estim_group(pattern, which='raw')
        estim_group(pattern, which='median_5')
        #estim_group(pattern, which='radon_5')
        
