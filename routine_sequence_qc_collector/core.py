import collections
import csv
import glob
import json
import logging
import os
import re
import shutil

from typing import Iterator, Optional

import routine_sequence_qc_collector.parsers as parsers
import routine_sequence_qc_collector.samplesheet as samplesheet


def create_output_dirs(config):
    """
    Create output directories if they don't exist.

    :param config: Application config.
    :type config: dict[str, object]
    :return: None
    :rtype: None
    """
    base_outdir = config['output_dir']
    output_dirs = [
        base_outdir,
        os.path.join(base_outdir, 'multiqc'),
        os.path.join(base_outdir, 'fastqc'),
        os.path.join(base_outdir, 'library-qc'),
        os.path.join(base_outdir, 'species-abundance'),
        os.path.join(base_outdir, 'bracken-species-abundances'),
    ]
    for output_dir in output_dirs:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)    


def find_latest_routine_sequence_qc_output(analysis_dir):
    """
    Find the latest routine sequence QC output directory, for a given run's analysis directory.

    :param analysis_dir: Analysis directory.
    :type analysis_dir: str
    :return: Path to latest routine sequence QC output directory.
    :rtype: str
    """
    routine_sequence_qc_output_dir_glob = "routine-sequence-qc-v*-output"
    routine_sequence_qc_output_dirs = glob.glob(os.path.join(analysis_dir, routine_sequence_qc_output_dir_glob))
    latest_routine_sequence_qc_output_dir = None
    if len(routine_sequence_qc_output_dirs) > 0:
        latest_routine_sequence_qc_output_dir = os.path.abspath(routine_sequence_qc_output_dirs[-1])

    return latest_routine_sequence_qc_output_dir


def find_analysis_dirs(config, check_complete=True):
    """
    Find all analysis directories.

    :param config: Application config.
    :type config: dict[str, object]
    :param check_complete: Check if analysis is complete.
    :type check_complete: bool
    :return: Analysis directory.
    :rtype: Iterator[Optional[dict[str, str]]]
    """
    miseq_run_id_regex = "\d{6}_M\d{5}_\d+_\d{9}-[A-Z0-9]{5}"
    nextseq_run_id_regex = "\d{6}_VH\d{5}_\d+_[A-Z0-9]{9}"
    analysis_by_run_dir = config['analysis_by_run_dir']
    subdirs = os.scandir(analysis_by_run_dir)
    
    for subdir in subdirs:
        run_id = subdir.name
        matches_miseq_regex = re.match(miseq_run_id_regex, run_id)
        matches_nextseq_regex = re.match(nextseq_run_id_regex, run_id)
        sequencer_type = None
        if matches_miseq_regex:
            sequencer_type = 'miseq'
        elif matches_nextseq_regex:
            sequencer_type = 'nextseq'
        not_excluded = run_id not in config['excluded_runs']
        ready_to_collect = False
        if check_complete:
            latest_routine_sequence_qc_output = find_latest_routine_sequence_qc_output(subdir)
            if latest_routine_sequence_qc_output is not None and os.path.exists(latest_routine_sequence_qc_output):
                routine_sequence_qc_analysis_complete = os.path.exists(os.path.join(latest_routine_sequence_qc_output, 'pipeline_complete.json'))
                ready_to_collect = routine_sequence_qc_analysis_complete
        else:
            ready_to_collect = True

        conditions_checked = {
            "is_directory": subdir.is_dir(),
            "matches_illumina_run_id_format": ((matches_miseq_regex is not None) or (matches_nextseq_regex is not None)),
            "not_excluded": not_excluded,
            "ready_to_collect": ready_to_collect,
        }
        conditions_met = list(conditions_checked.values())

        analysis_directory_path = os.path.abspath(subdir.path)
        analysis_dir = {
            "path": analysis_directory_path,
            "sequencer_type": sequencer_type,
        }
        if all(conditions_met):
            logging.info(json.dumps({
                "event_type": "analysis_directory_found",
                "sequencing_run_id": run_id,
                "analysis_directory_path": analysis_directory_path
            }))

            yield analysis_dir
        else:
            logging.debug(json.dumps({
                "event_type": "directory_skipped",
                "analysis_directory_path": os.path.abspath(subdir.path),
                "conditions_checked": conditions_checked
            }))
            yield None

            
def find_runs(config):
    """
    Finda all runs that have routine sequence QC data.

    :param config: Application config.
    :type config: dict[str, object]
    :return: List of runs. Keys: ['run_id', 'sequencer_type']
    :rtype: list[dict[str, str]]
    """
    logging.info(json.dumps({"event_type": "find_runs_start"}))
    runs = []
    all_analysis_dirs = sorted(list(os.listdir(config['analysis_by_run_dir'])))
    all_run_ids = filter(lambda x: re.match('\d{6}_[VM]', x) != None, all_analysis_dirs)
    for run_id in all_run_ids:
        if run_id in config['excluded_runs']:
            continue

        sequencer_type = None
        if re.match('\d{6}_M\d{5}_', run_id):
            sequencer_type = 'miseq'
        elif re.match('\d{6}_VH\d{5}_', run_id):
            sequencer_type = 'nextseq'

        analysis_dir = os.path.join(config['analysis_by_run_dir'], run_id)
        latest_routine_sequence_qc_output_dir = find_latest_routine_sequence_qc_output(analysis_dir)
        qc_check_complete_file = os.path.join(latest_routine_sequence_qc_output_dir, 'qc_check_complete.json')

        qc_check_info = {}
        if os.path.exists(qc_check_complete_file):
            with open(qc_check_complete_file, 'r') as f:
                qc_check_info = json.load(f)

        check_metrics = qc_check_info.get('checked_metrics', [])
        if os.path.exists(os.path.join(latest_routine_sequence_qc_output_dir, 'pipeline_complete.json')):
            run = {
                'run_id': run_id,
                'sequencer_type': sequencer_type,
                'run_qc_check': {
                    'checked_metrics': check_metrics,
                    'overall_qc_pass_fail': qc_check_info.get('overall_pass_fail', None),
                }
            }
            runs.append(run)

    logging.info(json.dumps({
        "event_type": "find_runs_complete"
    }))

    return runs


def scan(config: dict[str, object]) -> Iterator[Optional[dict[str, str]]]:
    """
    Scanning involves looking for all existing runs and...

    :param config: Application config.
    :type config: dict[str, object]
    :return: A run directory to analyze, or None
    :rtype: Iterator[Optional[dict[str, object]]]
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    for analysis_dir in find_analysis_dirs(config):    
        yield analysis_dir


def infer_species(config: dict[str, object], species_abundance, project_id):
    """
    :return: Name of inferred species.
    :rtype: str
    """
    species_abundance_keys = [
        'abundance_5_',
        'abundance_4_',
        'abundance_3_',
        'abundance_2_',
        'abundance_1_',
    ]
    inferred_species = None
    if 'projects' in config and project_id in config['projects']:
        project = config['projects'][project_id]
        if 'fixed_genome_size' in project and project['fixed_genome_size']:
            inferred_species = project.get('project_species_name', None)
        else:
            greatest_fraction_total_reads = 0.0
            for k in species_abundance_keys:
                if species_abundance[k + 'name'] != 'Homo sapiens' and species_abundance[k + 'fraction_total_reads'] > greatest_fraction_total_reads:
                    inferred_species = species_abundance[k + 'name']
                    greatest_fraction_total_reads = species_abundance[k + 'fraction_total_reads']
    else:
        greatest_fraction_total_reads = 0.0
        for k in species_abundance_keys:
            if (k + 'name') in species_abundance and species_abundance[k + 'name'] != 'Homo sapiens' and species_abundance[k + 'fraction_total_reads'] > greatest_fraction_total_reads:
                inferred_species = species_abundance[k + 'name']
                greatest_fraction_total_reads = species_abundance[k + 'fraction_total_reads']

    return inferred_species


def get_percent_reads_by_species_name(species_abundance, species_name):
    """
    """
    percent_reads = None
    species_abundance_keys = [
        'abundance_1_',
        'abundance_2_',
        'abundance_3_',
        'abundance_4_',
        'abundance_5_',
    ]

    for k in species_abundance_keys:
        try:
            if species_abundance[k + 'name'] == species_name:
                percent_reads = 100 * species_abundance[k + 'fraction_total_reads']
        except KeyError as e:
            pass

    return percent_reads        
    
    
def collect_outputs(config: dict[str, object], analysis_dir: Optional[dict[str, str]]):
    """
    Collect all routine sequence QC outputs for a specific analysis dir.

    :param config: Application config.
    :type config: dict[str, object]
    :param analysis_dir: Analysis dir. Keys: ['path', 'sequencer_type']
    :type analysis_dir: dict[str, str]
    :return: None
    :rtype: None
    """
    run_id = os.path.basename(analysis_dir['path'])
    logging.info(json.dumps({"event_type": "collect_outputs_start", "sequencing_run_id": run_id, "analysis_dir_path": analysis_dir['path']}))

    latest_routine_sequence_qc_output_path = find_latest_routine_sequence_qc_output(analysis_dir['path'])

    parsed_samplesheet_src_file = os.path.join(latest_routine_sequence_qc_output_path, 'parse_sample_sheet', 'sample_sheet.json')
    # If we can't find the parsed SampleSheet then we don't have a
    # Simple way to get Sample IDs and Project IDs. 
    if not os.path.exists(parsed_samplesheet_src_file):
        logging.error(json.dumps({'event_type': 'find_parsed_samplesheet_failed', 'sequencing_run_id': run_id, 'parsed_samplesheet_path': parsed_samplesheet_src_file}))
        return None

    libraries_by_library_id = {}
    with open(parsed_samplesheet_src_file, 'r') as f:
        samplesheet = json.load(f)
        if analysis_dir['sequencer_type'] == 'miseq':
            samplesheet_key = 'data'
        elif analysis_dir['sequencer_type'] == 'nextseq':
            samplesheet_key = 'cloud_data'
        else:
            logging.error(json.dumps({'event_type': 'find_parsed_samplesheet_failed', 'sequencing_run_id': run_id, 'parsed_samplesheet_path': parsed_samplesheet_src_file}))
            return None

        for sample in samplesheet[samplesheet_key]:
            if 'sample_id' in sample and re.match("S\d{1,4}$", sample['sample_id']):
                if 'sample_name' in sample:
                    library_id = sample['sample_name'].strip().replace('_', '-').replace('.', '-')
            else:
                library_id = sample['sample_id'].strip().replace('_', '-').replace('.', '-')
            if 'project_name' in sample:
                samplesheet_project_id = sample['project_name']
            elif 'sample_project' in sample:
                samplesheet_project_id = sample['sample_project']
            else:
                samplesheet_project_id = ""

            library = {
                'library_id': library_id,
                'samplesheet_project_id': samplesheet_project_id,
            }

            if samplesheet_project_id in config['projects']:
                project = config['projects'][samplesheet_project_id]
                if 'translated_project_id' in project and project['translated_project_id'] != '':
                    translated_project_id = project['translated_project_id']
                    library['translated_project_id'] = translated_project_id

            if 'translated_project_id' in library:
                library['project_id'] = library['translated_project_id']
            else:
                library['project_id'] = library['samplesheet_project_id']

            # This gets pretty verbose, even for debugging.
            # Un-comment during development if needed
            # logging.debug(json.dumps({
            #     'event_type': 'library_parsed_from_samplesheet',
            #     'parsed_samplesheet_path': parsed_samplesheet_src_file,
            #     'sequencing_run_id': run_id,
            #     'library_id': library_id,
            #     'project_id': project_id,
            # }))

            libraries_by_library_id[library_id] = library

    # species-abundance
    species_abundance_by_library_id = {library_id: {'library_id': library_id, 'project_id': libraries_by_library_id[library_id]['project_id']} for library_id in libraries_by_library_id.keys()}
    species_abundance_dst_file = os.path.join(config['output_dir'], "species-abundance", run_id + "_species_abundance.json")
    if not os.path.exists(species_abundance_dst_file):
        species_abundance_src_file = os.path.join(latest_routine_sequence_qc_output_path, 'abundance_top_n', 'top_5_abundances_species.csv')
        if os.path.exists(species_abundance_src_file):
            with open(species_abundance_src_file, 'r') as f:
                reader = csv.DictReader(f, dialect='unix')
                for row in reader:
                    library_id = row['sample_id']
                    if library_id in species_abundance_by_library_id:
                        for n in range(1, 6):
                            species_name = None
                            species_name_key = 'abundance_' + str(n) + '_name'
                            fraction_total_reads = None
                            fraction_total_reads_key = 'abundance_' + str(n) + '_fraction_total_reads'
                            species_name = row.get(species_name_key, None)
                            species_abundance_by_library_id[library_id][species_name_key] = species_name
                            try:
                                fraction_total_reads = float(row.get(fraction_total_reads_key, None))
                            except ValueError as e:
                                logging.error(json.dumps({'event_type': 'collect_species_abundance_metric_failed', 'metric': fraction_total_reads_key, 'sequencing_run_id': run_id, 'library_id': library_id}))
                            species_abundance_by_library_id[library_id][fraction_total_reads_key] = fraction_total_reads

        with open(species_abundance_dst_file, 'w') as f:
            json.dump(list(species_abundance_by_library_id.values()), f, indent=2)

        logging.info(json.dumps({
            "event_type": "write_species_abundance_complete",
            "run_id": run_id,
            "dst_file": species_abundance_dst_file
        }))

    if not os.path.exists(os.path.join(config['output_dir'], "bracken-species-abundances", run_id)):
        os.makedirs(os.path.join(config['output_dir'], "bracken-species-abundances", run_id))
        
    for library_id in libraries_by_library_id.keys():
        bracken_abundances_dst_file = os.path.join(config['output_dir'], "bracken-species-abundances", run_id, library_id + "_bracken_species_abundances.tsv")
        if not os.path.exists(bracken_abundances_dst_file):
            bracken_abundances_src_file = os.path.join(latest_routine_sequence_qc_output_path, 'bracken', library_id + '_Species_bracken_abundances_adjusted.tsv')
            if os.path.exists(bracken_abundances_src_file):
                shutil.copyfile(bracken_abundances_src_file, bracken_abundances_dst_file)
                logging.debug(json.dumps({
                    "event_type": "copy_bracken_abundances_complete",
                    "run_id": run_id,
                    "src_file": bracken_abundances_src_file,
                    "dst_file": bracken_abundances_dst_file
                }))

    # library-qc
    library_qc_dst_file = os.path.join(config['output_dir'], "library-qc", run_id + "_library_qc.json")
    if not os.path.exists(library_qc_dst_file):
        basic_qc_stats_src_file = os.path.join(latest_routine_sequence_qc_output_path, 'basic_qc_stats', 'basic_qc_stats.csv')
        if os.path.exists(basic_qc_stats_src_file):
            with open(basic_qc_stats_src_file, 'r') as f:
                reader = csv.DictReader(f, dialect='unix')
                for row in reader:
                    library_id = row['sample_id']
                    if library_id in libraries_by_library_id:
                        project_id = libraries_by_library_id[library_id]['project_id']
                        total_bases = None
                        percent_bases_above_q30 = None
                        inferred_species = None
                        inferred_species = infer_species(config, species_abundance_by_library_id[library_id], project_id)
                        if inferred_species is not None:
                            logging.debug(json.dumps({'event_type': 'library_species_inferred', 'sequencing_run_id': run_id, 'library_id': library_id, 'inferred_species': inferred_species}))
                            libraries_by_library_id[library_id]['inferred_species_name'] = inferred_species
                            percent_inferred_species = get_percent_reads_by_species_name(species_abundance_by_library_id[library_id], inferred_species)
                            if percent_inferred_species is None:
                                logging.error(json.dumps({"event_type": "collect_library_qc_metric_failed", "metric": "inferred_species_percent", 'library_id': library_id, 'inferred_species': inferred_species}))
                            libraries_by_library_id[library_id]['inferred_species_percent'] = percent_inferred_species
                            if 'known_species' in config and inferred_species in config['known_species']:
                                inferred_species_genome_size = config['known_species'][inferred_species]['genome_size_mb']
                                libraries_by_library_id[library_id]['inferred_species_genome_size_mb'] = inferred_species_genome_size
                        else:
                            logging.debug(json.dumps({'event_type': 'library_species_inference_failed', 'sequencing_run_id': run_id, 'library_id': library_id, 'inferred_species': inferred_species}))
                        try:
                            total_bases = int(row.get('total_bases', None))
                        except ValueError as e:
                            logging.error(json.dumps({'event_type': 'collect_library_qc_metric_failed', 'metric': 'total_bases', 'sequencing_run_id': run_id, 'library_id': library_id}))
                        libraries_by_library_id[library_id]['total_bases'] = total_bases
                        if all(k in libraries_by_library_id[library_id] for k in ['total_bases', 'inferred_species_genome_size_mb', 'inferred_species_percent']):
                            total_bases = libraries_by_library_id[library_id]['total_bases']
                            genome_size = libraries_by_library_id[library_id]['inferred_species_genome_size_mb'] * 1000000
                            species_percent = libraries_by_library_id[library_id]['inferred_species_percent']
                            if all([total_bases, genome_size, species_percent]):
                                libraries_by_library_id[library_id]['inferred_species_estimated_depth'] = (total_bases * (species_percent / 100)) / genome_size
                        try:
                            percent_bases_above_q30 = float(row.get('percent_bases_above_q30', None))
                        except ValueError as e:
                            logging.error(json.dumps({'event_type': 'collect_library_qc_metric_failed', 'metric': 'percent_bases_above_q30', 'sequencing_run_id': run_id, 'library_id': library_id}))
                        libraries_by_library_id[library_id]['percent_bases_above_q30'] = percent_bases_above_q30

        with open(library_qc_dst_file, 'w') as f:
            json.dump(list(libraries_by_library_id.values()), f, indent=2)

        logging.info(json.dumps({
            "event_type": "write_library_qc_complete",
            "run_id": run_id,
            "dst_file": library_qc_dst_file
        }))

    if not os.path.exists(os.path.join(config['output_dir'], "fastqc", run_id)):
        os.makedirs(os.path.join(config['output_dir'], "fastqc", run_id))

    for library_id in libraries_by_library_id.keys():
        # fastqc
        for read_type in ['R1', 'R2']:
            fastqc_src_file = os.path.join(latest_routine_sequence_qc_output_path, 'fastqc', '_'.join([library_id, read_type, 'fastqc']), 'fastqc_report.html')
            fastqc_dst_file = os.path.join(config['output_dir'], "fastqc", run_id, '_'.join([library_id, read_type, 'fastqc.html']))
            if not os.path.exists(fastqc_src_file):
                logging.warn(json.dumps({
                    "event_type": "copy_fastqc_failed",
                    "run_id": run_id,
                    "src_file": fastqc_src_file,
                    "dst_file": fastqc_dst_file,
                }))
            if os.path.exists(fastqc_src_file) and not os.path.exists(fastqc_dst_file):
                shutil.copyfile(fastqc_src_file, fastqc_dst_file)
                logging.debug(json.dumps({
                    "event_type": "copy_fastqc_complete",
                    "run_id": run_id,
                    "src_file": fastqc_src_file,
                    "dst_file": fastqc_dst_file
                })) 

    # multiqc
    multiqc_src_file = os.path.join(latest_routine_sequence_qc_output_path, 'multiqc', 'multiqc_report.html')
    multiqc_dst_file = os.path.join(config['output_dir'], "multiqc", run_id + "_multiqc.html")
    if not os.path.exists(multiqc_src_file):
        logging.warn(json.dumps({
            "event_type": "copy_multiqc_failed",
            "run_id": run_id,
            "src_file": multiqc_src_file,
            "dst_file": multiqc_dst_file
        }))
    if os.path.exists(multiqc_src_file) and not os.path.exists(multiqc_dst_file):
        shutil.copyfile(multiqc_src_file, multiqc_dst_file)
        logging.info(json.dumps({
            "event_type": "copy_multiqc_complete",
            "run_id": run_id,
            "src_file": multiqc_src_file,
            "dst_file": multiqc_dst_file
        }))


    logging.info(json.dumps({"event_type": "collect_outputs_complete", "sequencing_run_id": run_id, "analysis_dir_path": analysis_dir['path']}))
