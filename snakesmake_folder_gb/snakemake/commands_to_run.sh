#!/bin/bash

conda env create --name snakemake-chip_pipeline --file envs/chip_pipeline.yaml

conda activate snakemake-chip_pipeline

cd chip_pipeline

snakeamke --cores 4 --job 4

# snakemake activate command (need to make rules)