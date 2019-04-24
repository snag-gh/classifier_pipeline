from pprint import pprint as pp
from collections import defaultdict
import os
import pandas as pd
import datetime
from time import strftime

PIPELINE = os.environ['PIPELINE_HOME']
RUNDIR = os.environ['RUN_DIR']
OUTDIR = os.environ['OUT_DIR']
currentdate = os.environ['DATE']
configfile: PIPELINE + '/config/cluster.json'
GROUP = config['group']
VERSION = config['version']

SAMPLES = {}
TARGETS = []

samplesheet = RUNDIR + '/' + config['RUNID'] + '/Sample_Sheet.csv'
df0 = pd.read_csv(samplesheet, delimiter = ',', skiprows = 7, skip_blank_lines = True, dtype = {'Sample_Name': object, 'Sample_Well': object, 'Sample_Plate': object, 'Sample_Group': object, 'Pool_ID': object, 'Sentrix_ID': object, 'Sentrix_Position': object, 'Material_Type': object, 'Gender': object, 'Surgical_Case': object})
df = df0.dropna(axis=0, how='all')
df = pd.DataFrame(df)
for i, j, k in zip(df.Sample_Name.values, df.Sentrix_ID.values, df.Sentrix_Position.values):
    SAMPLES[i] = j + '_' + k
    

#pp(df)                           
#pp(SAMPLES)

TARGETS.append(OUTDIR + '/QCreports/' + config['RUNID'] + '/' + config['RUNID'] + '.qcReport.pdf')


for sample, sentrix_id in SAMPLES.items():
    TARGETS.append(OUTDIR + '/ClassifierReports/' + sample + '_' + sentrix_id + '/' + sample + '_Report_' + sentrix_id + '_run_' + currentdate + '.pdf')

#pp(TARGETS)

RUN = config['RUNID']
onstart:
    print('Started pipeline')
    shell("ssh biowulf.nih.gov \"echo 'Classifier pipeline {VERSION} started on run: {RUN}' | mutt -s 'Classifier Pipeline: {RUN}' `whoami`\@mail.nih.gov\"")

onsuccess:
    print('Pipeline completed')
    shell("find {OUTDIR}/ClassifierReports {OUTDIR}/QCreports/{RUN} -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {OUTDIR}/ClassifierReports {OUTDIR}/QCreports/{RUN} \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("ssh biowulf.nih.gov \"echo 'Classifier pipeline {VERSION} completed successfully on run: {RUN}' | mutt -s 'Classifier Pipeline: {RUN}' `whoami`\@mail.nih.gov\"")
    shell("find .snakemake/ logs {RUN}.yaml \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find .snakemake/ logs {RUN}.yaml -group $USER -exec chgrp -f {GROUP} {{}} \;")

onerror:
    print('An error occured')
    shell("find {OUTDIR}/ClassifierReports {OUTDIR}/QCreports/{RUN} -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {OUTDIR}/ClassifierReports {OUTDIR}/QCreports/{RUN} \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("ssh biowulf.nih.gov \"echo 'Classifier pipeline {VERSION} has errored on run: {RUN}' | mutt -s 'Classifier Pipeline: {RUN}' `whoami`\@mail.nih.gov\"")
    shell("find .snakemake/ logs {RUN}.yaml \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find .snakemake/ logs {RUN}.yaml -group $USER -exec chgrp -f {GROUP} {{}} \;")

rule all:
    input:
        TARGETS

rule QC:
    input:
        sampleSheet = RUNDIR + '/{runid}/Sample_Sheet.csv'

    output:
        main = OUTDIR + '/QCreports/{runid}/{runid}.qcReport.pdf',
        suppl = OUTDIR + '/QCreports/{runid}/{runid}.supplementary_plots.pdf'

    params:
        inputPath = RUNDIR + '/{runid}',
        rulename = 'QC.{runid}',
        resources = config['QC']

    shell:
        '''
        module load R/3.5.2
        {PIPELINE}/scripts/methylArrayQC.R {params.inputPath} {input.sampleSheet} {output.main} {output.suppl}
        '''

rule classifier:
    input:
        sampleSheet = RUNDIR + '/{runid}/Sample_Sheet.csv'

    output:
        temp(OUTDIR + '/ClassifierReports/{sample}_{runid}_{position}/NCIreport.Rmd'),
        pdf = OUTDIR + '/ClassifierReports/{sample}_{runid}_{position}/{sample}_Report_{runid}_{position}_run_' + currentdate + '.pdf'

    params:
        outdir = OUTDIR + '/ClassifierReports/{sample}_{runid}_{position}/',
        rulename = 'classifier.{runid}.{position}.{sample}',
        resources = config['classifier']

    shell:
        '''
        module load R/3.5
        module load pandoc/2.1.1
        module load tex/2018
        cp {PIPELINE}/scripts/classifier/NCIreport.Rmd {params.outdir}
        {PIPELINE}/scripts/classifier.R {input.sampleSheet} {output.pdf} {wildcards.sample} {wildcards.runid}_{wildcards.position} {PIPELINE} {VERSION}
        '''  
