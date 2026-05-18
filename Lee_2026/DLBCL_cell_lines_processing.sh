#!/bin/bash
#
# Request one node:
#SBATCH --array=1-16
#
# Specify one task:
#SBATCH --ntasks=1
#
# Number of processors for single task needed for use case (example):
#SBATCH --cpus-per-task=8
#
# Wall clock limit:
#SBATCH --time=05:00:00
#
## Command(s) to run (example):
export IGDATA=~/test-igblast
curr_folder=DLBCL_cell_lines/output
presto_output_folder=presto_assembled_first_output
igblast_output_folder=igblast_output_assembled_first

j_primers="primers/AbSeq/AbSeq_R1_Mouse_IG_Primers.fasta"
v_primers="primers/AbSeq/AbSeq_R2_TS.fasta"


while getopts "m" OPT; do
    case "$OPT" in
    m)  miseq=true
        ;;
    \?) echo -e "Invalid option: -${OPTARG}" >&2
        exit 1
        ;;
    :)  echo -e "Option -${OPTARG} requires an argument" >&2
        exit 1
        ;;
    esac
done

if [ "$miseq" = true ]; then
    echo "Running MiSeq Presto pipeline"
    j_primers="primers/MiSeq/CPrimers.fasta"
    v_primers="primers/MiSeq/VPrimers.fasta"
    outdir="DLBCL_cell_lines/MiSeq_output"
else
    echo "Running NEB Presto pipeline"
    outdir="DLBCL_cell_lines/NEB_output"
fi

presto_output_folder=${outdir}_presto_output
igblast_output_folder=${outdir}_igblast_output

doigblast() {
    local name=$1
    ~/changeo/bin/AssignGenes.py igblast -s ${curr_folder}/${presto_output_folder}/${name}.fastq -b ~/test-igblast    --organism mouse --loci ig --format blast --outdir ${curr_folder}/${igblast_output_folder}
    ~/changeo/bin/MakeDb.py igblast -i ${curr_folder}/${igblast_output_folder}/${name}_igblast.fmt7 -s ${curr_folder}/${presto_output_folder}/${name}.fastq    -r ~/test-germlines/imgt/mouse/vdj/imgt_mouse_*.fasta --extended
}

dopresto() {
    local cell_line=$1
    local sample=$2

    if [ "$miseq" = true ]; then
        doMiseqPresto $cell_line $sample
    else
        doNEBpresto $cell_line $sample
    fi
}

doNEBpresto() {
    # Assign the arguments to variables
    local cell_line=$1
    local sample=$2

    # get read file names and make unique sample names
    local r1="${sample}_L001_R1_001"
    local r2="${sample}_L001_R2_001"

    # run through the presto pipeline
    sh ~/presto-abseq.sh -1 ${curr_folder}/${r1}.fastq -2 ${curr_folder}/${r2}.fastq -j ${j_primers} -v ${v_primers} -o ${curr_folder}/${presto_output_folder}/ -p 8 -b "~/bin/igblastn" -n ${cell_line} -r ~/test-igblast/fasta/imgt_mouse_ig_v.fasta || exit 1
}

doMiseqPresto() {
    mkdir -p ${curr_folder}/${presto_output_folder}/log
    cd ${curr_folder}/${presto_output_folder}
    local sample=$2
    local cell_line=$1

    # get read file names and make unique sample names
    local r2="${sample}_L001_R1_001"
    local r1="${sample}_L001_R2_001"

    AssemblePairs.py sequential -1 ${curr_folder}/${r1}.fastq -2 ${curr_folder}/${r2}.fastq --coord illumina --outname ${cell_line} --log ${curr_folder}/${presto_output_folder}/log/AP_${cell_line}.log --outdir . -r ${home}/test-igblast/fasta/imgt_mouse_ig_v.fasta --aligner blastn --exec ~/bin/igblastn --nproc 8 --rc tail || exit 1
    FilterSeq.py quality -s ${curr_folder}/${presto_output_folder}/${cell_line}_assemble-pass.fastq -q 20 --outname ${cell_line} --log ${curr_folder}/${presto_output_folder}/log/FS_${cell_line}.log --outdir . --nproc 8
    # masking primers
    MaskPrimers.py score -s ${curr_folder}/${presto_output_folder}/${cell_line}_quality-pass.fastq -p ${v_primers}  --mode mask --pf VPRIMER --outname ${cell_line}-FWD --log ${curr_folder}/${presto_output_folder}/log/MPV_${cell_line}.log --outdir . --nproc 8
    MaskPrimers.py score -s ${curr_folder}/${presto_output_folder}/${cell_line}-FWD_primers-pass.fastq -p ${j_primers}  --mode mask --pf CPRIMER --outname ${cell_line}-REV --log ${curr_folder}/${presto_output_folder}/log/MPC_${cell_line}.log --revpr --outdir . --nproc 8


    # collapse seqs, split by dup count, parse headers, parse log
    CollapseSeq.py -s ${curr_folder}/${presto_output_folder}/${cell_line}-REV_primers-pass.fastq -n 20 --inner --outname ${cell_line} --outdir ${curr_folder}/${presto_output_folder} --log ${curr_folder}/${presto_output_folder}/log/CS_${cell_line}.log 
    SplitSeq.py group -s ${curr_folder}/${presto_output_folder}/${cell_line}_collapse-unique.fastq -f DUPCOUNT --num 2 --outname ${cell_line} --outdir ${curr_folder}/${presto_output_folder} 
    ParseHeaders.py table -s ${curr_folder}/${presto_output_folder}/${cell_line}_atleast-2.fastq -f ID DUPCOUNT --outname ${cell_line} --outdir ${curr_folder}/${presto_output_folder}
}

# List of cell lines to process
cell_lines=(
    "Bcl2_S1"
    "CBP1_S10"
    "CBP2_S11"
    "CMBP3_S12"
    "CMBP7_S13"
    "EBMP2_S8"
    "EBMP3_S9"
    "EZBP1_S4"
    "EZBP6_S5"
    "MCDP2_S14"
    "MCDP4_S15"
    "MCDP6_S16"
    "SN2F_S7"
    "SN2M_S6"
    "WP2_S2"
    "WP6_S3"
)

n=${SLURM_ARRAY_TASK_ID}-1


cell_line_sample=${cell_lines[$n]}

# get the first part of the string "cell_line" before the underscore
cell_line=$(echo $cell_line_sample | cut -d'_' -f1)

echo "Processing cell line: $cell_line"
dopresto $cell_line $cell_line_sample
presto_name=${cell_line}-final_collapse-unique_atleast-2
doigblast ${presto_name}
igblast_name=${presto_name}_igblast_db-pass
