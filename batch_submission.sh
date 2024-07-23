#!/bin/bash

for f in data/config_files/*.yaml; do
	echo "$f"
	cp -rf $f data/config.yaml
<<<<<<< HEAD
	snakemake --use-conda --cores all -R create_input
=======
	snakemake --use-conda --cores all
>>>>>>> v.1/v.1
done
