#!/usr/bin/env bash
# Script to run the pRESTO initial assembly and annotation on COVID data
#
# Author:  Ishita Singh, modified from script by Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    04/01/2024
#
# Arguments:
#   -1  Read 1 FASTQ sequence file (sequence beginning with the C-region or J-segment).
#   -2  Read 2 FASTQ sequence file (sequence beginning with the leader or V-segment).
#   -j  Read 1 FASTA primer sequences (C-region or J-segment).
#       Defaults to forward_primer.fa
#   -v  Read 2 FASTA primer sequences (template switch or V-segment).
#       Defaults to reverse_primers.fa.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to a directory matching the sample identifier in the current working directory.
#   -x  The mate-pair coordinate format of the raw data.
#       Defaults to illumina.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -1  Read 1 FASTQ sequence file.\n" \
            "     Sequence beginning with the C-region or J-segment)."
    echo -e "  -2  Read 2 FASTQ sequence file.\n" \
            "     Sequence beginning with the leader or V-segment)."
    echo -e "  -j  Read 1 FASTA primer sequences.\n" \
            "    Defaults to forward_primer.fa"
    echo -e "  -v  Read 2 FASTA primer or template switch sequences.\n" \
            "     Defaults to reverse_primers.fa"
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the read 1 filename."
    echo -e "  -o  Output directory. Will be created if it does not exist.\n" \
            "     Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -x  The mate-pair coordinate format of the raw data.\n" \
            "     Defaults to illumina."
    echo -e "  -h  This message."
}

# Argument validation variables
R1_READS_SET=false
R2_READS_SET=false
FWD_PRIMERS_SET=false
REV_PRIMERS_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
COORD_SET=false

# Get commandline arguments
while getopts "1:2:j:v:r:f:n:o:x:p:h" OPT; do
    case "$OPT" in
    1)  R1_READS=$OPTARG
        R1_READS_SET=true
        ;;
    2)  R2_READS=$OPTARG
        R2_READS_SET=true
        ;;
    j)  FWD_PRIMERS=$OPTARG
        FWD_PRIMERS_SET=true
        ;;
    v)  REV_PRIMERS=$OPTARG
        REV_PRIMERS_SET=true
        ;;
    n)  OUTNAME=$OPTARG
        OUTNAME_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    x)  COORD=$OPTARG
        COORD_SET=true
        ;;
    h)  print_usage
        exit
        ;;
    \?) echo -e "Invalid option: -${OPTARG}" >&2
        exit 1
        ;;
    :)  echo -e "Option -${OPTARG} requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if required arguments are not provided
if ! ${R1_READS_SET} || ! ${R2_READS_SET}; then
    echo -e "You must specify both read files using the -1 and -2 options." >&2
    exit 1
fi

# Set unspecified arguments
if ! ${OUTNAME_SET}; then
    OUTNAME=$(basename ${R1_READS} | sed 's/\([0-9]\{3\}[A-Z]\).*/\1/')
fi

if ! ${OUTDIR_SET}; then
    OUTDIR=${OUTNAME}
fi

if ! ${COORD_SET}; then
    COORD="illumina"
fi

# Check R1 reads
if [ -e ${R1_READS} ]; then
    R1_READS=$(realpath ${R1_READS})
else
    echo -e "File ${R1_READS} not found." >&2
    exit 1
fi

# Check R2 reads
if [ -e ${R2_READS} ]; then
    R2_READS=$(realpath ${R2_READS})
else
    echo -e "File ${R2_READS} not found." >&2
    exit 1
fi

# Check R1 primers
if ! ${FWD_PRIMERS_SET}; then
    FWD_PRIMERS="../forward_primer.fa"
elif [ -e ${FWD_PRIMERS} ]; then
    FWD_PRIMERS=$(realpath ${FWD_PRIMERS})
else
    echo -e "File ${FWD_PRIMERS} not found." >&2
    exit 1
fi

# Check R2 primers
if ! ${REV_PRIMERS_SET}; then
    REV_PRIMERS="../reverse_primers.fa"
elif [ -e ${REV_PRIMERS} ]; then
    REV_PRIMERS=$(realpath ${REV_PRIMERS})
else
    echo -e "File ${REV_PRIMERS} not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
FILTER_LOWQUAL=true

# FilterSeq run parameter
FS_QUAL=20

# MaskPrimers run parameters
MP_MAXLEN=1000
MP_MAXERROR=0.1


# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="${OUTNAME}_logs"
TABDIR="${OUTNAME}_parse_files"
PIPELINE_LOG="${LOGDIR}/pipeline-assemble.log"
ERROR_LOG="${LOGDIR}/pipeline-assemble.err"
mkdir -p ${LOGDIR}
mkdir -p ${TABDIR}
echo '' > $PIPELINE_LOG
echo '' > $ERROR_LOG

# Check for errors
check_error() {
    if [ -s $ERROR_LOG ]; then
        echo -e "ERROR:"
        cat $ERROR_LOG | sed 's/^/    /'
        exit 1
    fi
}

# Start
PRESTO_VERSION=$(python3 -c "import presto; print('%s-%s' % (presto.__version__, presto.__date__))")
echo -e "IDENTIFIER: ${OUTNAME}"
echo -e "DIRECTORY: ${OUTDIR}"
echo -e "PRESTO VERSION: ${PRESTO_VERSION}"
echo -e "\nSTART"
STEP=0

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
AssemblePairs.py align -1 ${R2_READS} -2 ${R1_READS}  \
    --coord $COORD --rc tail --log "${LOGDIR}/assemble.log" --outname "${OUTNAME}" --outdir . \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Remove low quality reads
if $FILTER_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
    #OUTPREFIX="$(printf '%02d' $STEP)--${OUTNAME}"
    FilterSeq.py quality -s "${OUTNAME}_assemble-pass.fastq" -q $FS_QUAL \
        --outname "${OUTNAME}" --log "${LOGDIR}/quality.log" \
        >> $PIPELINE_LOG  2> $ERROR_LOG
    MP_FILE="${OUTNAME}_quality-pass.fastq"
    check_error
else
    MP_FILE="${OUTNAME}_assemble-pass.fastq"
fi

# Masking reads (align and then score)
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s $MP_FILE -p $FWD_PRIMERS --mode mask \
    --maxlen $MP_MAXLEN --maxerror $MP_MAXERROR --pf VPRIMER \
    --log "${LOGDIR}/fwd-primer.log" --outname "${OUTNAME}_fwd" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s "${OUTNAME}_fwd_primers-pass.fastq" -p $REV_PRIMERS --mode cut \
    --start 0 --maxerror $MP_MAXERROR --pf CPRIMER --revpr --outname "${OUTNAME}_rev" \
    --log "${LOGDIR}/rev-primer.log" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Note: the next couple of steps are not logged in the workflow
# Removal of duplicate sequencies and filtering to repeated sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s "${OUTNAME}_rev_primers-pass.fastq" -n 20 --inner --uf CPRIMER \
    --cf VPRIMER --act set --outname "${OUTNAME}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Filtering to repeated sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq"
SplitSeq.py group -s "${OUTNAME}_collapse-unique.fastq" -f DUPCOUNT --num 2 --outname "${OUTNAME}" \
  >> $PIPELINE_LOG 2> $ERROR_LOG



# Annotation table
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders"
ParseHeaders.py table -s "${OUTNAME}_atleast-2.fastq" -f ID DUPCOUNT CPRIMER VPRIMER --outname "${OUTNAME}" --outdir ${TABDIR}\
>> $PIPELINE_LOG 2> $ERROR_LOG


# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l "${LOGDIR}/assemble.log" \
    -f ID LENGTH OVERLAP ERROR PVALUE \
    --outdir ${TABDIR} > /dev/null 2> $ERROR_LOG &
if $FILTER_LOWQUAL; then
    ParseLog.py -l "${LOGDIR}/quality.log" -f ID QUALITY --outdir ${TABDIR} \
         > /dev/null 2> $ERROR_LOG &
fi
ParseLog.py -l "${LOGDIR}/fwd-primer.log" "${LOGDIR}/rev-primer.log" -f ID PRIMER ERROR \
    --outdir ${TABDIR} > /dev/null 2> $ERROR_LOG &
wait
check_error


# Zip or delete intermediate and log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
FILTER_FILES="$(basename ${OUTDIR}_atleast-2.fastq)\|$(basename ${FWD_PRIMERS})\|$(basename ${REV_PRIMERS})"
FILTER_FILES+="\|Amplicon-sequencing\|Amplicon-sequencing.zip"
TEMP_FILES=$(ls *.fastq | grep -v ${FILTER_FILES})
if $ZIP_FILES; then
    tar -zcf "${OUTNAME}_log_files.tar.gz" "${LOGDIR}"
    tar -zcf "${OUTNAME}_parse_files.zip" "${TABDIR}"
    tar -zcf "${OUTNAME}_temp_files.tar.gz" $TEMP_FILES
fi
if $DELETE_FILES; then
    rm $TEMP_FILES
    rm -r "${LOGDIR}"
    rm -r "${TABDIR}"
fi


# End
printf "DONE\n\n"
cd ../

