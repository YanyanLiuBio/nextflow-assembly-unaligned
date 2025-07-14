#!/bin/bash

path=s3://seqwell-analysis/20250604_MiSeq-i100/plasmid_assembly_per_base/02Jun2025_prep_01_DR008a_96_well_format

/software/nextflow-align/nextflow run \
main.nf \
--path $path \
-resume -bg 
