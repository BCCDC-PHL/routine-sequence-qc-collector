import collections
import re
import csv


def parse_artic_qc(artic_qc_path, run_id):
    """
    """
    output = []
    all_input_fields = [
        "sample_name",
        "pct_N_bases",
        "pct_covered_bases",
        "longest_no_N_run",
        "num_aligned_reads",
        "fasta",
        "bam",
        "qc_pass",
    ]

    int_fields = [
        'longest_no_N_run',
        'num_aligned_reads',
    ]

    float_fields = [
        'pct_N_bases',
        'pct_covered_bases',
    ]

    boolean_fields = [
        'qc_pass',
    ]
    
    with open(artic_qc_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            qc = collections.OrderedDict()
            for field in all_input_fields:
                if row[field] == 'NA':
                    qc[field] = None
                elif field == 'sample_name':
                    library_id = row[field]
                    qc['library_id'] = library_id
                    plate_id = None
                    if library_id.startswith('POS') or library_id.startswith('NEG'):
                        try:
                            plate_id = re.search("\d+", library_id.split('-')[2]).group(0)
                        except AttributeError as e:
                            pass
                    else:
                        try:
                            plate_id = re.search("\d+", library_id.split('-')[1]).group(0)
                        except AttributeError as e:
                            pass
                    if plate_id:
                        qc['plate_id'] = int(plate_id)
                    qc['run_id'] = run_id
                elif field == 'pct_N_bases':
                    pass
                elif field == 'qc_pass':
                    pass
                elif field == 'pct_covered_bases':
                    try:
                        qc['genome_completeness'] = float(row[field])
                    except ValueError as e:
                        qc[field] = None
                elif field in int_fields:
                    try:
                        qc[field] = int(row[field])
                    except ValueError as e:
                        qc[field] = None
                elif field in float_fields:
                    try:
                        qc[field] = float(row[field])
                    except ValueError as e:
                        qc[field] = None
                elif field in boolean_fields:
                    try:
                        qc[field] = bool(row[field])
                    except ValueError as e:
                        qc[field] = None
                elif field == 'run_name':
                    qc['plate_id'] = int(row[field].split('_')[-1])
                    qc['run_id'] = '_'.join(row[field].split('_')[0:-1])
                elif field == 'qc_pass':
                    qc[field] = row[field].split(',')
                else:
                    qc[field] = row[field]
            output.append(qc)

    return output


def parse_ncov_tools_summary_qc(ncov_tools_summary_qc_path):
    """
    """
    output = []
    all_input_fields = [
        "sample",
        "run_name",
        "num_consensus_snvs",
        "num_consensus_n",
        "num_consensus_iupac",
        "num_variants_snvs",
        "num_variants_indel",
        "num_variants_indel_triplet",
        "mean_sequencing_depth",
        "median_sequencing_depth",
        "qpcr_ct",
        "collection_date",
        "num_weeks",
        "scaled_variants_snvs",
        "genome_completeness",
        "qc_pass",
        "lineage",
        "lineage_notes",
        "watch_mutations",
    ]

    int_fields = [
        'num_consensus_snvs',
        'num_consensus_n',
        'num_consensus_iupac',
        'num_variants_snvs',
        'num_variants_indel',
        'num_variants_indel_triplet',
        'median_sequencing_depth',
        'num_weeks',
    ]

    float_fields = [
        'mean_sequencing_depth',
        'qpcr_ct',
        'scaled_variants_snvs',
        'genome_completeness',
    ]
    with open(ncov_tools_summary_qc_path, 'r') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        for row in reader:
            qc = collections.OrderedDict()
            for field in all_input_fields:
                if row[field] == 'NA':
                    qc[field] = None
                elif field in int_fields:
                    try:
                        qc[field] = int(row[field])
                    except ValueError as e:
                        qc[field] = None
                elif field in float_fields:
                    try:
                        qc[field] = float(row[field])
                    except ValueError as e:
                        qc[field] = None
                elif field == 'sample':
                    qc['library_id'] = row[field]
                elif field == 'run_name':
                    qc['plate_id'] = int(row[field].split('_')[-1])
                    qc['run_id'] = '_'.join(row[field].split('_')[0:-1])
                elif field == 'qc_pass':
                    qc[field] = row[field].split(',')
                else:
                    qc[field] = row[field]
            output.append(qc)

    return output


def parse_amplicon_depth_bed(amplicon_depth_bed_path):
    """
    """
    output = []
    all_input_fields = [
        'reference_name',
        'start',
        'end',
        'amplicon_id',
        'pool',
        'strand',
        'mean_depth',
    ]

    int_fields = [
        'start',
        'end',
        'pool',
    ]

    float_fields = [
        'mean_depth',
    ]

    with open(amplicon_depth_bed_path, 'r') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        for row in reader:
            amplicon = collections.OrderedDict()
            for field in all_input_fields:
                if row[field] == 'NA':
                    amplicon[field] = None
                elif field in int_fields:
                    try:
                        amplicon[field] = int(row[field])
                    except ValueError as e:
                        amplicon[field] = None
                elif field in float_fields:
                    try:
                        amplicon[field] = float(row[field])
                    except ValueError as e:
                        amplicon[field] = None
                elif field == 'reference_name':
                    pass
                elif field == 'amplicon_id':
                    try:
                        amplicon['amplicon_num'] = int(row[field].split('_')[-1])
                    except ValueError as e:
                        amplicon['amplicon_num'] = None
                else:
                    amplicon[field] = row[field]
            output.append(amplicon)

    return output
