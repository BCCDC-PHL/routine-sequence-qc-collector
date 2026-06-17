import glob
import json
import logging
import os
import re

import routine_sequence_qc_collector.instrument as instrument

def find_samplesheet_for_run(run_id, sequencer_output_dirs):
    """
    """
    samplesheets = []
    samplesheet_path = None
    instrument_type = instrument.determine_instrument_type(run_id)
    logging.debug(json.dumps({"event_type": "determined_instrument_type", "run_id": run_id, "instrument_type": instrument_type}))

    if instrument_type == "unknown":
        return samplesheet_path

    sequencer_run_dir = ""
    if instrument_type == 'miseq':
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

    elif instrument_type == 'nextseq':
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

    elif instrument_type == 'i100':
        for sequencer_output_dir in sequencer_output_dirs:
            if not re.search('i100', sequencer_output_dir):
                continue
            for run_dir in os.listdir(sequencer_output_dir):
                run_dir = os.path.abspath(os.path.join(sequencer_output_dir, run_dir))
                if os.path.basename(run_dir) == run_id:
                    sequencer_run_dir = run_dir
                    logging.debug(json.dumps({"event_type": "found_sequencer_output_dir", "run_id": run_id, "sequencer_output_dir": sequencer_run_dir}))
                    exit()
            

    if len(samplesheets) == 1:
        samplesheet_path = samplesheets[0]

    return samplesheet_path
