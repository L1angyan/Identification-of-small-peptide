bowtie2="~/software/miniconda3/envs/common/bin/bowtie2"
bowtie2_index="~/project/peptide/ref/bowtie2_index/genome"
STAR="~/software/miniconda3/envs/common/bin/STAR"
STAR_index="~/project/peptide/ref/STAR_index"
samtools="~/command/software/samtools"
sambamba="~/software/miniconda3/envs/common/bin/sambamba"
bamCoverage="~/software/miniconda3/envs/common/bin/bamCoverage"
RiboCode_onestep="~/software/miniconda3/envs/ribocode/bin/RiboCode_onestep"
GTF="~/project/peptide/ref/Zea_mays.B73_RefGen_v4.50.chr.gtf"
fasta="~/project/peptide/ref/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
genomesize="2200000000"

for i in $(cat zzsample.txt|sort|uniq);do
# merge replicate bam
$sambamba sort -t 8 mapping/${i}.rep1.Aligned.toTranscriptome.out.bam
$sambamba sort -t 8 mapping/${i}.rep2.Aligned.toTranscriptome.out.bam
$sambamba merge -t 8 mapping/$i.merge.bam $(ls mapping/${i}.rep?.Aligned.toTranscriptome.out.sorted.bam | tr "\n" " ")
$sambamba sort -t 8 mapping/$i.merge.bam

# RiboCode pipeline
$RiboCode_onestep \
    --gtf $GTF --fasta $fasta --rpf_mapping_file mapping/$i.merge.sorted.bam \
    --longest-orf no --min-AA-length 0 --pval-cutoff 0.01 \
    --output-name res_merge/RiboCode_ORFs_res.$i --output-gtf
done
