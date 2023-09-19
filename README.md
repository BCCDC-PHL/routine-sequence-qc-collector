# routine-sequence-qc-collector

## Usage

```
usage: routine-sequence-qc-collector [-h] [-c CONFIG] [--log-level LOG_LEVEL]

options:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
  --log-level LOG_LEVEL
```

## Configuration

A `config-template.json` file is provided in this repo. The tool expects a json-formatted config file with these fields:

```json
{
    "analysis_by_run_dir": "/path/to/routine-sequence-qc/analysis_by_run",
    "excluded_runs_list": "/path/to/excluded_runs.csv",
    "projects_definition_file": "/path/to/projects.csv",
    "known_species_list": "/path/to/known_species.csv",
    "scan_interval_seconds": 3600,
    "output_dir": "/path/to/routine-sequence-qc-collector/data"
}
```