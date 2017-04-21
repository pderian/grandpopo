"""Runs the GPP data pre-processing steps.

- frame extraction from source video;
- image rectification and filtering;
- output as image files along with json metadata.

Written by P. DERIAN, 2016-2017.
www.pierrederian.net
"""
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

if __name__=="__main__":

    action = 'decoder'

    if action=='decoder':
        # the area centered over sensors
        pattern = '20140313/1[2-6]/*.mp4' #sensors
        process_group(pattern, params=PARAMS_COMP60)
        #process_group(pattern, params=PARAMS_COMP30)
        # the wide field for flash rip monitoring (TGRS paper)
        #pattern = '20140312/12/*.mp4' #flash rip
        #process_group(pattern, params=PARAMS_RIP120)
        # the wide field for swash (COASTAL DYN paper)
        #pattern = '20140313/15/*.mp4' # [TMP]
        #process_group(pattern, params=PARAMS_SWASH125)

#### Note on igrida runs:
# oarsub -I -p "host='igrida-abacus.irisa.fr' AND gpu='YES'" -l walltime=00:30:00
#
# module load cuda
# /udd/pderian/resources/TyphoonPack_v2/build_yupana/CuTyphoon/cutyphoon_server -h -nD 6 -nE 6 -nM 7 -nPyr 1 -nI 3 -r 1 -a 1e-3 -port 40000 -dev 1 > /temp_dd/igrida-fs1/pderian/SCRATCH/gpp/estim/timeseries_60m/typhoon_log.txt 2>&1 &
# python ~/scripts/Python/grandpopo/estimator.py 40000 -o gpp/estim/timeseries_60m/ --probe --no-fields -c 'timeseries_60m estim with nD 6 nE 6 nM 7 nPyr 1 nI 3 r 1 a 1e-3' -f gpp/plage/timeseries_60m/20140313/20140313_12/20140313_1201/*.json
