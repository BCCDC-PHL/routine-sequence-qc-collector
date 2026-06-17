import re

MISEQ_RUN_ID_REGEX = "\\d{6}_M\\d{5}_\\d+_\\d{9}-[A-Z0-9]{5}"
NEXTSEQ_RUN_ID_REGEX = "\\d{6}_VH\\d{5}_\\d+_[A-Z0-9]{9}"
I100_RUN_ID_REGEX = "\\d{8}_SH\\d{5}_\\d+_[A-Z0-9]{10}-[A-Z0-9]{3}"

RUN_ID_REGEX_BY_INSTRUMENT_TYPE = {
    'nextseq': NEXTSEQ_RUN_ID_REGEX,
    'miseq': MISEQ_RUN_ID_REGEX,
    'i100': I100_RUN_ID_REGEX,
}


def matches_a_valid_run_id_regex(run_id: str):
    """
    Determine if the run ID matches any valid Run ID regex.
    """
    matches = False
    for regex_instrument_type, regex in RUN_ID_REGEX_BY_INSTRUMENT_TYPE.items():
        if re.match(regex, run_id):
            matches = True

    return matches


def determine_instrument_type(run_id: str):
    """
    Determine the instrument type based on matching the run ID aginst a regex.
    :param run_id: The run ID
    :type run_id: str
    :return: The instrument type: one of: 'miseq', 'nextseq', 'i100', 'unknown'
    :rtype: str
    """
    instrument_type = 'unknown'

    for regex_instrument_type, regex in RUN_ID_REGEX_BY_INSTRUMENT_TYPE.items():
        matches_regex = re.match(regex, run_id)
        if matches_regex:
            instrument_type = regex_instrument_type

    return instrument_type


