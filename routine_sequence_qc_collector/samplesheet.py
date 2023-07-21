import glob
import json
import logging
import os
import re


def find_samplesheet_for_run(run_id, sequencer_output_dirs):
    """
    """
    samplesheets = []
    samplesheet_path = None
    sequencer_type = None
    if re.match('\d{6}_M', run_id):
        sequencer_type = 'miseq'
    elif re.match('\d{6}_VH', run_id):
        sequencer_type = 'nextseq'
    else:
        return samplesheet_path

    sequencer_run_dir = ""
    if sequencer_type == 'miseq':
        for sequencer_output_dir in sequencer_output_dirs:
            if not re.search('miseq', sequencer_output_dir):
                continue
            for run_dir in os.listdir(sequencer_output_dir):
                if os.path.basename(run_dir) == run_id:
                    sequencer_run_dir = os.path.join(sequencer_output_dir, run_dir)

        if os.path.exists(os.path.join(sequencer_run_dir, 'Alignment_1')):
            # Run is 'new-style' MiSeq Output directory
            demultiplexing_output_dirs = os.listdir(os.path.join(sequencer_run_dir, 'Alignment_1'))
            most_recent_demultiplexing_outdir = sorted(demultiplexing_output_dirs)[-1]
            samplesheets = [os.path.join(most_recent_demultiplexing_outdir, 'SampleSheetUsed.csv')]
        else:
            # Run is 'old-style' MiSeq Output directory
            standard_samplesheet_path = os.path.join(sequencer_run_dir, 'SampleSheet.csv')
            if os.path.exists(standard_samplesheet_path):
                samplesheets = [standard_samplesheet_path]
                logging.debug(json.dumps({"event_type": "found_samplesheets", "run_id": run_id, "sequencer_run_dir": sequencer_run_dir, "samplesheet_paths": samplesheets}))
            else:
                samplesheets = glob.glob(os.path.join(sequencer_run_dir, 'SampleSheet*.csv'))
                logging.debug(json.dumps({"event_type": "found_samplesheets", "run_id": run_id, "sequencer_run_dir": sequencer_run_dir, "samplesheet_paths": samplesheets}))
            
    elif sequencer_type == 'nextseq':
        logging.debug(json.dumps({"event_type": "determined_sequencer_type", "run_id": run_id, "sequencer_type": sequencer_type}))
        for sequencer_output_dir in sequencer_output_dirs:
            if not re.search('nextseq', sequencer_output_dir):
                continue
            for run_dir in os.listdir(sequencer_output_dir):
                run_dir = os.path.abspath(os.path.join(sequencer_output_dir, run_dir))
                if os.path.basename(run_dir) == run_id:
                    sequencer_run_dir = run_dir
                    logging.debug(json.dumps({"event_type": "found_sequencer_output_dir", "run_id": run_id, "sequencer_output_dir": sequencer_run_dir}))

        if os.path.exists(os.path.join(sequencer_run_dir, 'Analysis')):
            demultiplexing_output_dirs = os.listdir(os.path.join(sequencer_run_dir, 'Analysis'))
            most_recent_demultiplexing_outdir = os.path.join(sequencer_run_dir, 'Analysis', sorted(demultiplexing_output_dirs)[-1])
            logging.debug(json.dumps({"event_type": "determined_most_recent_demultiplexing_outdir", "run_id": run_id, "demultiplexing_outdir": most_recent_demultiplexing_outdir}))
            samplesheets = glob.glob(os.path.join(most_recent_demultiplexing_outdir, 'Data', 'SampleSheet*.csv'))

    if len(samplesheets) == 1:
        samplesheet_path = samplesheets[0]

    return samplesheet_path


def count_covid19_production_samples_in_samplesheet(samplesheet_path, sequencer_type):
    """
    """
    num_covid19_production_samples = 0
    if sequencer_type == 'nextseq':
        with open(samplesheet_path, 'r') as f:
            for row in f:
                if re.search(',covid-19_production,', row):
                    num_covid19_production_samples += 1
    elif sequencer_type == 'miseq':
        covid_library_id_regexes = [
            'S\d{1,3},[ER]\d{10},',                                  # Container ID only
            'S\d{1,3},[FHSTW]\d{6},',                                # Foreign Container ID
            'S\d{1,3},X\d{5},',                                      # Foreign Container ID
            'S\d{1,3},[ER]\d{10}-\d{1,4}-[A-Z0-9]{1,2}-[A-H]\d{2},', # Container ID-Plate Num-Index Set ID-Well
            'S\d{1,3},POS',                                          # Positive control
            'S\d{1,3},NEG',                                          # Negative control
        ]
        with open(samplesheet_path, 'r') as f:
            for row in f:
                library_id_matches = []
                for r in covid_library_id_regexes:
                    if re.match(r, row):
                        library_id_matches.append(re.match(r, row).group(0))
                    else:
                        library_id_matches.append(None)
                if any(library_id_matches):
                    num_covid19_production_samples += 1

    return num_covid19_production_samples
