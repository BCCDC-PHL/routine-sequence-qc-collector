#!/usr/bin/env python

import argparse
import datetime
import json
import logging
import os
import time

import routine_sequence_qc_collector.config
import routine_sequence_qc_collector.core as core

from routine_sequence_qc_collector.logging_config import configure_logging

DEFAULT_SCAN_INTERVAL_SECONDS = 3600.0

log = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config')
    parser.add_argument('--log-level')
    args = parser.parse_args()

    configure_logging(args.log_level)

    config = {}
    scan_interval = DEFAULT_SCAN_INTERVAL_SECONDS

    quit_when_safe = False

    while(True):
        if quit_when_safe:
            exit(0)
        try:
            if args.config:
                try:
                    config = routine_sequence_qc_collector.config.load_config(args.config)
                    # Un-comment below to debug config parsing
                    # config['excluded_runs'] = list(config['excluded_runs'])
                    # print(json.dumps(config, indent=2))
                    # exit()
                    log.info({"event_type": "config_loaded", "config_file": os.path.abspath(args.config)})
                except json.decoder.JSONDecodeError as e:
                    # If we fail to load the config file, we continue on with the
                    # last valid config that was loaded.
                    log.error({"event_type": "load_config_failed", "config_file": os.path.abspath(args.config)})

            core.create_output_dirs(config)

            scan_start_timestamp = datetime.datetime.now()

            runs = core.find_runs(config)
            runs_output_file = os.path.join(config['output_dir'], 'runs.json')
            with open(runs_output_file, 'w') as f:
                json.dump(runs, f, indent=2)
            log.info({"event_type": "write_runs_file_complete", "runs_file": runs_output_file})

            for run in core.scan(config):
                if run is not None:
                    try:
                        config = routine_sequence_qc_collector.config.load_config(args.config)
                        log.info({"event_type": "config_loaded", "config_file": os.path.abspath(args.config)})
                    except json.decoder.JSONDecodeError as e:
                        log.error({"event_type": "load_config_failed", "config_file": os.path.abspath(args.config)})
                    core.collect_outputs(config, run)
                if quit_when_safe:
                    exit(0)
            scan_complete_timestamp = datetime.datetime.now()
            scan_duration_delta = scan_complete_timestamp - scan_start_timestamp
            scan_duration_seconds = scan_duration_delta.total_seconds()
            scan_interval_seconds = config.get('scan_interval_seconds', None)
            if scan_interval_seconds:
                next_scan_timestamp = datetime.datetime.now() + datetime.timedelta(seconds=scan_interval_seconds)
            log.info({"event_type": "scan_complete", "scan_duration_seconds": scan_duration_seconds, "timestamp_next_scan": str(next_scan_timestamp.isoformat())})

            if quit_when_safe:
                exit(0)

            if "scan_interval_seconds" in config:
                try:
                    config['scan_interval_seconds'] = float(str(config['scan_interval_seconds']))
                except ValueError as e:
                    config['scan_interval_seconds'] = DEFAULT_SCAN_INTERVAL_SECONDS
            else:
                    config['scan_interval_seconds'] = DEFAULT_SCAN_INTERVAL_SECONDS
            time.sleep(config['scan_interval_seconds'])
        except KeyboardInterrupt as e:
            log.info({"event_type": "quit_when_safe_enabled"})
            quit_when_safe = True

if __name__ == '__main__':
    main()
