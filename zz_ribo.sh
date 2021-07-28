cd /public/home/yliang/ref/rat

# GTF文件升级
GTFupdate Rattus_norvegicus.Rnor_6.0.104.gtf > Rattus_updated.gtf

# 准备RiboCode所有文件
prepare_transcripts -g Rattus_updated.gtf \
-f Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
-o RiboCode_annot

# star建索引
mkdir star_index
module load STAR/2.7.3a

bsub -J index -n 10 -R span[hosts=1] -o %J.out -e %J.err -q normal \
"STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ./star_index \
--genomeFastaFiles ./Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
--sjdbGTFfile ./Rattus_updated.gtf"

# tRNA,rRNA,snoRNA序列建索引，以便去除Ribo-seq中的这些序列
mkdir ribo_filter
awk '!/^>/{printf "%s", $0; n="\n"}/^>/{print n $0; n=""}\
END{printf "%s", n}' $1|grep -A1 -E "tRNA|rRNA|snoRNA" \
> ./ribo_filter/nc_filter.fa
bsub -J filter -n 5 -R span[hosts=1] -o %J.out -e %J.err -q q2680v2 \
"bowtie2-build --threads 5 ./ribo_filter/nc_filter.fa ./ribo_filter/filter"

# 下载解压
for i in `cat zz_rat.download`;
do
bsub -J download -n 1 -R span[hosts=1] -o %J.out -e %J.err -q q2680v2 \
"prefetch $i -O ./rat; \
fastq-dump --split-3 --gzip $2/${i}.sra -O ./rat;"
done

# 去除adapter
for i in *fastq.gz;do echo -e $i,`dnapi.py $i`; done > adapter.txt
for i in `cat adapter.txt`;
do 
file=${i%,*}; 
out=${file%.fastq*}.fq.gz; 
adapter=${i#*,} 
bsub -J trim -n 5 -R span[hosts=1] -o %J.out -e %J.err -q q2680v2 \
"cutadapt -a $adapter -j 5 -m 20 -M 35 -o $out $file"; 
done

# mapping & identify ORFs from Ribo-seq data
for i in *.fq.gz;
do
echo ${i%.fq.gz}
done > zz_sample.txt
module load STAR/2.7.3a
for i in `cat zz_sample.txt`;
do
mkdir $i
cd $i
bsub -J ribo -n 10 -R span[hosts=1] -o %J.out -e %J.err -q high \
"bowtie2 -x ~/ref/ecoli/ribo_filter/filter -p 10 \
-U ../${i}.fq.gz -S ${i}.sam --un-gz ${i}_clean.fq.gz; \
rm ${i}.sam; \
STAR --outFilterType BySJout --runThreadN 10 --outFilterMismatchNmax 2 \
--genomeDir ~/ref/ecoli/star_index/ --readFilesCommand zcat \
--readFilesIn ${i}_clean.fq.gz --outFileNamePrefix ${i}_ --limitBAMsortRAM 1500000000 \
--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd; \
metaplots -a ~/ref/ecoli/RiboCode_annot -r ${i}_Aligned.toTranscriptome.out.bam; \
RiboCode -a ~/ref/ecoli/RiboCode_annot -c metaplots_pre_config.txt -l no -m 2 \
-p 0.01 -g -o ${i}_result;"
cd ../
done

# 整合所有结果
for i in `cat zz_sample.txt`;
do
mv ${i}/${i}_result.gtf ${i}/${i}_result.gtf ./zz_result
done

cd ./zz_result
for i in *result.txt;do
name=${i%.txt};
gtf=${name}.gtf;
cat $i >> z_all.txt;
cat $gtf >> zzz_all.gtf;
done
head -n1 z_all.txt > zzz_all.txt;
grep -v "ORF_ID" z_all.txt >> zzz_all.txt;
rm z_all.txt;

python3 zzz_merge.py
# 合并gtf
awk 'BEGIN{FS="\t"}{print ">rat"$30"\n"$29}' zzz_merge.txt > zz_rat_aa.fa
# 提取氨基酸序列
python3 zzz_trans_gtf.py
sed s/ORF/transcript/ zzz_peptide_ok.gtf > zzz_peptide_rat.gtf
# 提取核酸序列
gffread -w zzz_peptide_rat.fa \
-g ~/ref/rat/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
zzz_peptide_rat.gtf