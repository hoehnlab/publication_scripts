# this script was run in Docker container, using this terminal command:
# docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.5.0 bash

studies=("003O", "003W", "009O", "009W", "025O", "025W", "049O", "049W", "346O", "346W", "418O", "418W")
threads=4

for study in ${studies[@]}
do
    ofile="/data/for_igb/${study}_atleast-2.fastq"
    ofile2="/data/for_igb/${study}_atleast-2_reheader.fasta"
    ofile3="/data/for_igb/${study}_atleast-2_reheader_igblast.fmt7"
    echo $ofile
    echo $ofile2
    echo $ofile3

    ParseHeaders.py add -s $ofile -f sample -u $study --fasta

    # run IgBlast to align to V and J genes
    AssignGenes.py igblast -s $ofile2 \
        -b ./ --organism 'human' --loci ig --format blast \
        --nproc $threads
    
    MakeDb.py igblast -i $ofile3 -s $ofile2 \
        -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta \
         /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta \
        /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta \
         /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGKV.fasta \
         /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGKJ.fasta \
         /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGLV.fasta \
         /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGLJ.fasta 
done     