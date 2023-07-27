import json
import csv


def get_excluded_runs(config):
    """
    """
    excluded_runs = set()
    with open(config['excluded_runs_list'], 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                run_id = line.strip()
                excluded_runs.add(run_id)

    return excluded_runs


def get_projects(config):
    """
    Read the 'projects_definition_file' to get details for projects.

    :param config: Application config
    :type config: dict[str, object]
    :return: Project details. Keys: ['samplesheet_project_id', 'translated_project_id', 'project_species_name', 'fixed_genome_size', 'genome_size_mb']
    :rtype: dict[str, object]
    """
    projects = {}
    with open(config['projects_definition_file'], 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            project = {}
            samplesheet_project_id = row['samplesheet_project_id']
            project['samplesheet_project_id'] = samplesheet_project_id
            translated_project_id = row['translated_project_id']
            if translated_project_id != '':
                project['translated_project_id'] = translated_project_id

            if row['project_species_name'] != '':
                project['project_species_name'] = row['project_species_name']

            if row['project_species_taxid'] != '':
                project['project_species_taxid'] = row['project_species_taxid']

            if row['fixed_genome_size'].lower() == 'true':
                project['fixed_genome_size'] = True
            else:
                project['fixed_genome_size'] = False

            if row['genome_size_mb'] != '':
                try:
                    genome_size_mb = float(row['genome_size_mb'])
                    project['genome_size_mb'] = genome_size_mb
                except ValueError as e:
                    pass

            projects[samplesheet_project_id] = project
            if translated_project_id != '':
                projects[translated_project_id] = project

    return projects


def get_known_species(config):
    """
    """
    known_species = {}
    with open(config['known_species_list'], 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            species = {}
            taxid = row['ncbi_taxonomy_id']
            species['ncbi_taxonomy_id'] = taxid
            
            if row['species_name'] != '':
                species['species_name'] = row['species_name']

            if row['genome_size_mb'] != '':
                try:
                    genome_size_mb = float(row['genome_size_mb'])
                    species['genome_size_mb'] = genome_size_mb
                except ValueError as e:
                    pass

            if row['gc_percent'] != '':
                try:
                    gc_percent = float(row['gc_percent'])
                    species['gc_percent'] = gc_percent
                except ValueError as e:
                    pass

            if row['refseq_assembly_accession'] != '':
                species['refseq_assembly_accession'] = row['refseq_assembly_accession']

            known_species[taxid] = species
            if species['species_name'] != "":
                known_species[species['species_name']] = species

    return known_species


def load_config(config_path: str) -> dict[str, object]:
    """
    """
    with open(config_path, 'r') as f:
        config = json.load(f)

    if 'excluded_runs_list' in config:
        excluded_runs = get_excluded_runs(config)
        config['excluded_runs'] = excluded_runs
    else:
        config['excluded_runs'] = set()

    if 'projects_definition_file' in config:
        projects = get_projects(config)
        config['projects'] = projects
    else:
        config['projects'] = {}

    if 'known_species_list' in config:
        known_species = get_known_species(config)
        config['known_species'] = known_species
    else:
        config['known_species'] = {}

    return config
