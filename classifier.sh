#! /bin/bash

function usage() {
    echo "USAGE: $0 --runid [RUNID] --dryrun --outdir [optional: /path/to/output/dir] --pipeline [optional: /path/to/pipeline/dir]"
}
function fail() {
    echo "$@"
    exit 1
}

if [ $# -eq 0 ]; then
    usage
    exit 0
fi

RUN_DIR=/data/Compass/iScan_raw
OUT_DIR=/data/Compass/Methylation
PIPELINE_HOME=/data/Compass/Tools/classifier_pipeline
dryrun=
runid=

while [ "$1" != "" ]; do
    case $1 in
        --runid )	shift
			runid=$1
			;;
        --dryrun )	dryrun=1
			;;
        --outdir )	shift
			OUT_DIR=$1
			;;
        --pipeline )	shift
			PIPELINE_HOME=$1
			;;
        -h | --help )	usage
			exit
    esac
    shift
done

export RUN_DIR
export OUT_DIR
export PIPELINE_HOME
export DATE=`date +'%m%d%Y_%H%M%S'`

echo RUN_DIR: $RUN_DIR
echo OUT_DIR: $OUT_DIR
echo PIPELINE_HOME: $PIPELINE_HOME
echo DATETIME: $DATE

YAML=${runid}.yaml
echo $YAML
echo "RUNID:" > $YAML
echo "    '$runid'" >> $YAML
cat $YAML

export YAML 

if [ ! -d "$RUN_DIR/$runid" ]; then
    fail "Could not find RUN $runid at the input path $RUN_DIR/"
fi

module load snakemake/5.4.0 &> /dev/null || fail "Could not load module snakemake/5.4.0"

if [ "$dryrun" == '1' ];then
    #Dryrun
    echo "Dryrun"
    snakemake -nrp --nolock -k -j 3000 -s $PIPELINE_HOME/classifier.snakefile --configfile $YAML -d `pwd` 
else
    echo "Executing pipeline on RUN: $runid"
    if [ ! -d "logs" ]; then
        mkdir logs
        chgrp -f Compass logs
        chmod g+rwx logs
    fi
    sbatch -e pipeline.%j.%x.e -o pipeline.%j.%x.o --job-name=mnp_pipeline.$runid.$DATE --mem=1G --partition=ccr,norm --time=04:00:00 --cpus-per-task=1 $PIPELINE_HOME/submit.sh	
fi
