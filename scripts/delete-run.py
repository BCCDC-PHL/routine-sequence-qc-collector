#!/usr/bin/env python3

import argparse
import glob
import os
import shutil

def main(args):
    """
    Delete the run data for a given run ID.
    """
    subdirs = [
        'bracken-species-abundances',
        'fastqc',
        'library-qc',
        'multiqc',
        'species-abundance',
    ]
    for subdir in subdirs:

        path_glob = os.path.join(args.data_dir, subdir, args.run_id + '*')
        for file_or_dir in glob.glob(path_glob):
            if os.path.isfile(file_or_dir):
                os.remove(file_or_dir)
            elif os.path.isdir(file_or_dir):
                shutil.rmtree(file_or_dir)
            else:
                pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run-id', required=True)
    parser.add_argument('-d', '--data-dir', default='/data/analysis/routine-sequence-qc/site-data')
    args = parser.parse_args()
    main(args)
