task mutation_validator_preprocess {

    #Inputs and constants defined here
    String PAIRID
    File MAF
    File? WEXTUMOR
    File? WEXNORMAL
    File? WEXTUMORBAI
    File? WEXNORMALBAI
    File? WGSTUMOR
    File? WGSNORMAL
    File? WGSTUMORBAI
    File? WGSNORMALBAI
    File? RNATUMOR
    File? RNATUMORBAI
    File? TARGTUMOR
    File? TARGNORMAL
    File? TARGTUMORBAI
    File? TARGNORMALBAI
    File? LPTUMOR
    File? LPNORMAL
    File? LPTUMORBAI
    File? LPNORMALBAI
    File? OTUMOR
    File? ONORMAL
    File? OTUMORBAI
    File? ONORMALBAI

    String output_disk_gb
    String boot_disk_gb = "1000"
    String ram_gb = "8"
    String cpu_cores = "2"
    String preemption
    command {
python_cmd="
import subprocess,os
def run(cmd):
    print('about to run')
    print(cmd)
    print('')
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################

# split MAF into INDELs and SNPs

run('python /opt/src/filter_tsv.py -i \"${MAF}\"  -f Variant_Type -v \"SNP|DNP|TNP|MNP\" -e .snp.maf')

run('python /opt/src/filter_tsv.py -i \"${MAF}\"  -f Variant_Type -v \"INS|DEL\" -e .indel.maf')

#run('mkdir -p softlinked ')

run('ls -latr ')

MAF1=\"${MAF}\"
MAFSNP=os.path.basename(MAF1)+'.snp.maf'
MAFINDEL=os.path.basename(MAF1)+'.indel.maf'
RNATYPE='hg19'
WEXTBAM=\"${WEXTUMOR}\"
WEXNBAM=\"${WEXNORMAL}\"
WGSTBAM=\"${WGSTUMOR}\"
WGSNBAM=\"${WGSNORMAL}\"
RNATBAM=\"${RNATUMOR}\"
TARGTBAM=\"${TARGTUMOR}\"
TARGNBAM=\"${TARGNORMAL}\"
LPTBAM=\"${LPTUMOR}\"
LPNBAM=\"${LPNORMAL}\"
OTBAM=\"${OTUMOR}\"
ONBAM=\"${ONORMAL}\"
CWD = os.getcwd() 

PAIRID = \"${PAIRID}\"
PREPROCESSED_FILE = '%s.pileup_preprocessing.txt'%PAIRID

if os.path.exists(WEXTBAM):
    run('ln -sT ' + WEXTBAM + ' WEXT.bam')
    run('ln -sT ' + \"${WEXTUMORBAI}\" + ' WEXT.bam.bai')
    WEXTBAM = 'WEXT.bam' #os.path.basename(WEXTBAM)
#    run('ls -latrh ' + WEXTBAM +'*')

else:
    WEXTBAM='None'

if os.path.exists(WEXNBAM):
    run('ln -sT ' + WEXNBAM + ' WEXN.bam')
    run('ln -sT ' + \"${WEXNORMALBAI}\" + ' WEXN.bam.bai' )
    WEXNBAM = 'WEXN.bam' #os.path.basename(WEXNBAM)
#    run('ls -latrh ' + WEXTBAM+'*')
else:
    WEXNBAM='None'

if os.path.exists(WGSTBAM):
    run('ln -sT ' + WGSTBAM + ' WGST.bam')
    run('ln -sT ' + \"${WGSTUMORBAI}\"  + '  WGST.bam.bai')
    WGSTBAM = CWD + '/WGST.bam' #os.path.basename(WGSTBAM)
#    run('ls -latrh ' + WGSTBAM+'*')
else:
    WGSTBAM='None'

if os.path.exists(WGSNBAM):
    run('ln -sT ' + WGSNBAM + ' WGSN.bam')
    run('ln -sT ' + \"${WGSNORMALBAI}\" + ' WGSN.bam.bai')
    WGSNBAM = CWD + '/WGSN.bam' #os.path.basename(WGSNBAM)
#    run('ls -latrh ' + WGSNBAM+'*')
else:
    WGSNBAM='None'

if os.path.exists(RNATBAM):
    run('ln -sT ' + RNATBAM + ' RNAT.bam')
    run('ln -sT ' + \"${RNATUMORBAI}\" + ' RNAT.bam.bai' )
    RNATBAM = 'RNAT.bam'  #os.path.basename(RNATBAM)
#    run('ls -latrh ' + RNATBAM+'*')
    run('samtools view -H ' + RNATBAM + ' | grep SN:chr1 > RNACHECK.txt')
    with open('RNACHECK.txt', 'r') as f:
        RNACHECK_first_line = f.readline()
    if len(RNACHECK_first_line)>0:
        print('file is aligned with chr in front')
        RNATYPE='hg19-chr'
else:
    RNATBAM='None'

if os.path.exists(TARGTBAM):
    run('ln -sT ' + TARGTBAM + ' TARGT.bam')
    run('ln -sT ' + \"${TARGTUMORBAI}\" + ' TARGT.bam.bai')
    TARGTBAM = 'TARGT.bam'  #os.path.basename(TARGTBAM)
#    run('ls -latrh ' + TARGTBAM+'*')
else:
    TARGTBAM='None'

if os.path.exists(TARGNBAM):
    run('ln -sT ' + TARGNBAM + ' TARGN.bam' )
    run('ln -sT ' + \"${TARGNORMALBAI}\" + '  TARGN.bam.bai' )
    TARGNBAM = 'TARGN.bam' #os.path.basename(TARGNBAM)
#    run('ls -latrh ' + TARGNBAM+'*')
else:
    TARGNBAM='None'

if False and os.path.exists(LPTBAM):
    run('ln -sT ' + LPTBAM + ' .')
    run('ln -sT ' + \"${LPTUMORBAI}\" + ' .' )
    LPTBAM = os.path.basename(LPTBAM)
#    run('ls -latrh ' + LPTBAM+'*')
else:
    LPTBAM='None'

if False and os.path.exists(LPNBAM):
    run('ln -sT ' + LPNBAM + ' .')
    run('ln -sT ' + \"${LPNORMALBAI}\" + ' .')
    LPNBAM = os.path.basename(LPNBAM)
#    run('ls -latrh ' + LPNBAM+'*')
else:
    LPNBAM='None'

if False and os.path.exists(OTBAM):
    run('ln -sT ' + OTBAM + ' .')
    run('ln -sT ' + \"${OTUMORBAI}\" + ' .')
    OTBAM = os.path.basename(OTBAM)
#    run('ls -latrh ' + OTBAM+'*')
else:
    OTBAM='None'

if False and os.path.exists(ONBAM):
    run('ln -sT ' + ONBAM + ' .')
    run('ln -sT ' + \"${ONORMALBAI}\" + ' .')
    ONBAM = os.path.basename(ONBAM)
#    run('ls -latrh ' + ONBAM+'*')
else:
    ONBAM='None'


cmd1=' --mafsnp '+ MAFSNP + ' --mafindel ' + MAFINDEL + ' --wextumor ' + WEXTBAM + ' --wexnormal ' + WEXNBAM + ' --wgstumor ' +  WGSTBAM + ' --wgsnormal ' + WGSNBAM
cmd2=' --rnatumor ' + RNATBAM + ' --targetedtumor ' + TARGTBAM + ' --targetednormal ' + TARGNBAM + ' --lowpasstumor ' + LPTBAM + ' --lowpassnormal ' + LPNBAM 
cmd3=' --othertumor ' + OTBAM + ' --othernormal ' + ONBAM + ' --out ' + \"${PAIRID}\" + ' --rnatype ' + RNATYPE
cmd='python /opt/src/mutation_validator_preprocess.py ' + cmd1 + cmd2 + cmd3
print cmd
run(cmd)

cmd='mkdir snp_mv && cd snp_mv && python /opt/src/algutil/firehose_module_adaptor/run_module.py --module_libdir /opt/src/fh_MutationValidator \
--mutation.validator.preprocessed.file ../%s \
--maf_file_to_annotate ../%s \
--discovery_type.wgs_wex_rna_targeted wex \
--pair_id  %s.snp \
--print_discovery_counts  true \
--normal_coverage_threshold  0 \
--job.spec.memory 2'%(PREPROCESSED_FILE,MAFSNP,PAIRID)

run(cmd)

cmd='mkdir indel_mv && cd indel_mv && python /opt/src/algutil/firehose_module_adaptor/run_module.py --module_libdir /opt/src/fh_MutationValidator \
--mutation.validator.preprocessed.file ../%s \
--maf_file_to_annotate ../%s \
--discovery_type.wgs_wex_rna_targeted wex \
--pair_id  %s.indel \
--print_discovery_counts  true \
--normal_coverage_threshold  0 \
--job.spec.memory 2'%(PREPROCESSED_FILE,MAFINDEL,PAIRID) 

run(cmd)


import time
time.sleep(999999999)


#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File pileup_preprocessing_txt="${PAIRID}.pileup_preprocessing.txt"
        File validated_snp_maf="snp_mv/${PAIRID}.snp.validated.maf"
        File validated_indel_maf="indel_mv/${PAIRID}.indel.validated.maf"
        File dstat_log="dstat.log"
    }

    runtime {
        docker : "docker.io/stewart/mutation_validator_preprocess:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemption}"
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow mutation_validator_preprocess_workflow {
    call mutation_validator_preprocess
}
