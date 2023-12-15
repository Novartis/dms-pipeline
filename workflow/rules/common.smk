import os


def get_samples():
    """Returns the list of samples (sample_name)"""
    return list(
        pd.read_csv(config["samples_config"], sep="\t", header=0)["sample_name"]
    )


def get_orig_fastq(wildcards, mate="fastq_name_R1"):
    """Returns the paths to the input fastq files"""
    sample_sheet = pd.read_csv(
        config["samples_config"], sep="\t", index_col=None, header=0
    )
    assert mate in (
        "fastq_name_R1",
        "fastq_name_R2",
    ), f"Not valid arg, mate={mate}. Should be one of 'fastq_name_R2' or 'fastq_name_R2'"
    sample_query = "sample_name == '%s'" % wildcards.s
    orig_dir = sample_sheet.query(sample_query).fastq_path.tolist()[0]
    orig = sample_sheet.query(sample_query)[mate].tolist()[0]
    return os.path.join(orig_dir, orig)
