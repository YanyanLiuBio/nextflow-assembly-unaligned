#!/bin/bash

aws s3 ls s3://seqwell-analysis/20250516_Admera_Health/spades_busco/QUAST/ --recursive \
  | grep "all_alignments_.*-contigs.tsv" \
  | awk '{print $4}' \
  | while read path; do
      aws s3 cp "s3://seqwell-analysis/$path" ./;
    done