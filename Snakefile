shell.executable("/bin/bash")

import pandas as pd
import os, multiprocessing
import yaml
from snakemake import utils
from snakemake.utils import min_version
from snakemake import logging

min_version("6.0")
shell.prefix("set -euo pipefail;")

config = yaml.load(open("data/config.yaml", "r+"), Loader=yaml.FullLoader)
dependencies = yaml.load(open("src/dependencies.yaml", "r+"), Loader=yaml.FullLoader)

snakefiles = "src/"
include: snakefiles + "rules.py"

rule all:
    input:
        plot_raw = "results/{protein_name}/plots/{protein_name}_rawIntensities.pdf".format(protein_name=config['protein_name']),
        plot_kinetics = "results/{protein_name}/plots/{protein_name}_finalKinetics.pdf".format(protein_name=config['protein_name']),
        final_kinetics = "results/{protein_name}/{protein_name}_finalKinetics.csv".format(protein_name=config['protein_name']),
        chromatogram = "results/{protein_name}/{protein_name}_TICs.done".format(protein_name=config['protein_name']),
        coveragevalues = "results/{protein_name}/{protein_name}_coveragevalues.RData".format(protein_name=config['protein_name']),
        residuemap = "results/{protein_name}/{protein_name}_residuemap.gif".format(protein_name=config['protein_name'])
