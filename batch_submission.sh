#!/bin/bash

for f in data/config_files/*.yaml; do
	echo "$f"
	cp -rf $f data/config.yaml
	snakemake --use-conda --cores all
done
