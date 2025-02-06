cutadapt="~/software/miniconda3/envs/common/bin/cutadapt"
bowtie2="~/software/miniconda3/envs/common/bin/bowtie2"
bowtie2_index="~/project/peptide/ref/bowtie2_index/genome"
STAR="~/software/miniconda3/envs/common/bin/STAR"
STAR_index="~/project/peptide/ref/STAR_index"
samtools="~/command/software/samtools"
bamCoverage="~/software/miniconda3/envs/common/bin/bamCoverage"
RiboCode_onestep="~/software/miniconda3/envs/ribocode/bin/RiboCode_onestep"
GTF="~/project/peptide/ref/Zea_mays.B73_RefGen_v4.50.chr.gtf"
fasta="~/project/peptide/ref/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
genomesize="2200000000"

mkdir mapping
mkdir QC
mkdir res

for i in $(cat zzsample.txt);do
# optional
# $cutadapt -a CTGTAGGCACCATCAAT -j 16 --minimum-length 20 --maximum-length 35 -o QC/$i.fastq.gz raw/$i.fastq.gz

# rRNA removal
$bowtie2 \
    -x ${bowtie2_index} -p 16 --norc --un QC/$i.fq.gz -U QC/$i.fastq.gz -S QC/$i.rRNA.sam

# align to Transcriptome
$STAR \
    --outFilterType BySJout --runThreadN 16 --outFilterMismatchNmax 2 --genomeDir ${STAR_index} \
    --readFilesIn QC/$i.fq.gz --outFileNamePrefix mapping/$i. --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 --alignEndsType EndToEnd

# coverage
$samtools index -@ 16 mapping/${i}.Aligned.sortedByCoord.out.bam
$bamCoverage -b mapping/${i}.Aligned.sortedByCoord.out.bam -o mapping/${i}.bw \
               -p 8 --binSize 50 --normalizeUsing RPKM --effectiveGenomeSize $genomesize

# RiboCode pipeline
$RiboCode_onestep \
    --gtf $GTF --fasta $fasta --rpf_mapping_file mapping/$i.Aligned.toTranscriptome.out.bam \
    --longest-orf no --min-AA-length 0 --pval-cutoff 0.01 --output-name res/RiboCode_ORFs_res.$i --output-gtf
done
