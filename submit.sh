#!/bin/bash

snakemake -rp --nolock -k -j 3000 -s $PIPELINE_HOME/classifier.snakefile --configfile $YAML -d `pwd` --jobscript $PIPELINE_HOME/scripts/jobscript.sh --jobname {params.rulename}.{jobid} --cluster "sbatch -o logs/{params.rulename}.%j.o -e logs/{params.rulename}.%j.e {params.resources}"
