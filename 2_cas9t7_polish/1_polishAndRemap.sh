set -e

cd work

CPUS=40
CPUS2=8
PILON="java -Xmx20G -jar pilon/pilon.jar"

curl "https://tritrypdb.org/common/downloads/release-50/LmexicanaMHOMGT2001U1103/fasta/data/TriTrypDB-50_LmexicanaMHOMGT2001U1103_Genome.fasta" -o genome.fasta
cp ../../1_remap_Fiebig_et_al/work/genome.gff genome.gff

#Fetch illumina reads
gunzip ../56_S4_R1_001.fastq.gz -c > F.fq
gunzip ../56_S4_R2_001.fastq.gz -c > R.fq

# #https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
perl rcorrector/run_rcorrector.pl -t 12 -1 F.fq -2 R.fq 

python TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F.cor.fq -2 R.cor.fq -s c9t7 
rm F.cor.fq
rm R.cor.fq

TrimGalore-0.6.0/trim_galore --path_to_cutadapt=./cutadapt_venv/bin/cutadapt --paired --retain_unpaired --phred33 --output_dir trimgalore_c9t7 --length 36 -q 5 --stringency 1 -e 0.1 unfixrm_F.cor.fq unfixrm_R.cor.fq 
rm unfixrm_F.cor.fq
rm unfixrm_R.cor.fq
mv trimgalore_c9t7/unfixrm_F.cor_val_1.fq F.cor.fq
mv trimgalore_c9t7/unfixrm_R.cor_val_2.fq R.cor.fq
mv trimgalore_c9t7/unfixrm_F.cor_unpaired_1.fq U.cor.fq
cat trimgalore_c9t7/unfixrm_R.cor_unpaired_2.fq >> U.cor.fq
rm trimgalore_c9t7/unfixrm_R.cor_unpaired_2.fq

#Polish for non-IUPAC output
./bwa/bwa index genome.fasta 
./bwa/bwa mem -t $CPUS genome.fasta F.cor.fq R.cor.fq | ./samtools/bin/samtools sort -o alignment.sort.bam -@ $CPUS2 -m 4G 
./samtools/bin/samtools index alignment.sort.bam -@ $CPUS2 
$PILON --genome genome.fasta --frags alignment.sort.bam --output genome --outdir pilonOut --changes --vcf --tracks --fix snps,indels || echo
cp pilonOut/genome.fasta genome.polish.fasta 
sed -i 's/_pilon//g' genome.polish.fasta
cat pilonOut/genome.changes | nodejs ../pilonChanges.js >> polish.changes.txt

#Remap GFF, create GFF for polished genome
cat pilonOut/genome.changes | nodejs ../pilonRemapGFF.js genome.gff > genome.polish.gff 
nodejs ../combineFastaGFF.js genome.polish.gff genome.polish.fasta > genome.polish.combined.gff 
nodejs ../cdsGFFtoFASTA.js genome.polish.combined.gff nucl > genome.polish.orfs.fasta 
nodejs ../cdsGFFtoFASTA.js genome.polish.combined.gff prot > genome.polish.prots.fasta 

#Generate reference GFF, for reference genome
nodejs ../combineFastaGFF.js genome.gff genome.fasta > genome.combined.gff 
nodejs ../cdsGFFtoFASTA.js genome.combined.gff nucl > genome.orfs.fasta 
nodejs ../cdsGFFtoFASTA.js genome.combined.gff prot > genome.prots.fasta 

rm *.bam
rm *.bai
rm *.fq
