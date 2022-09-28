set -e

cd work

CPUS=40
MINIMAP2="minimap2/minimap2"
MINIASM="miniasm/miniasm"
RACON="racon/build/bin/racon"
PILON="java -Xmx16G -jar pilon/pilon.jar"
SAMTOOLS="./samtools/bin/samtools"
PORECHOP="./venv/bin/porechop"
BWA="./bwa/bwa"

#Fetch nanopore reads
cp ../nanopore.reads1d2.fastq reads.nanopore.fastq
cat ../nanopore.reads1d.fastq >> reads.nanopore.fastq

#Chop adapters
$PORECHOP --threads $CPUS --input reads.nanopore.fastq --output reads.nanoporeTrimmed.fastq
cat reads.nanoporeTrimmed.fastq | nodejs ../fastaHisto.js > reads.nanoporeTrimmed.histo.txt

#Genome assembly
#Align reads to reads
$MINIMAP2 -t$CPUS -x ava-ont reads.nanoporeTrimmed.fastq reads.nanoporeTrimmed.fastq | gzip -1 > reads.paf.gz
#Determine consensus in overlaps
$RACON --threads $CPUS --include-unpolished reads.nanoporeTrimmed.fastq reads.paf.gz reads.nanoporeTrimmed.fastq > reads.nanoporeTrimmedConsensus.fasta
#Align consensus reads to consensus reads
$MINIMAP2 -t$CPUS -x ava-ont reads.nanoporeTrimmedConsensus.fasta reads.nanoporeTrimmedConsensus.fasta | gzip -1 > reads.paf.gz
#Generate assembly
$MINIASM -f reads.nanoporeTrimmedConsensus.fasta reads.paf.gz > assembly.gfa
rm reads.paf.gz
#Make FASTA from GFA
awk '/^S/{print ">"$2"\n"$3}' assembly.gfa > assembly.fasta
cat assembly.fasta | nodejs ../fastaLeng.js > assembly.lengths.txt

#Fetch illumina reads
gunzip -c ../../cas9t7_polish/56_S4_R1_001.fastq.gz > F.fq
gunzip -c ../../cas9t7_polish/56_S4_R2_001.fastq.gz > R.fq
cat F.fq | nodejs ../fastaHisto.js > reads.illuminaF.histo.txt
cat R.fq | nodejs ../fastaHisto.js > reads.illuminaR.histo.txt

#https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
perl rcorrector/run_rcorrector.pl -t 12 -1 F.fq -2 R.fq

python TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F.cor.fq -2 R.cor.fq -s c9t7
rm F.cor.fq
rm R.cor.fq

TrimGalore-0.6.0/trim_galore --path_to_cutadapt=./venv/bin/cutadapt --paired --retain_unpaired --phred33 --output_dir trimgalore_c9t7 --length 36 -q 5 --stringency 1 -e 0.1 unfixrm_F.cor.fq unfixrm_R.cor.fq
rm unfixrm_F.cor.fq
rm unfixrm_R.cor.fq
mv trimgalore_c9t7/unfixrm_F.cor_val_1.fq F.cor.fq
mv trimgalore_c9t7/unfixrm_R.cor_val_2.fq R.cor.fq
mv trimgalore_c9t7/unfixrm_F.cor_unpaired_1.fq U.cor.fq
cat trimgalore_c9t7/unfixrm_R.cor_unpaired_2.fq >> U.cor.fq
rm trimgalore_c9t7/unfixrm_R.cor_unpaired_2.fq

#Iterative pilon polish
ITER=0
rm polish.changes.txt || echo
echo "IUPAC polishing" > polish.changes.txt
cp assembly.fasta assembly.polished$ITER.fasta
mkdir pilonOut
ITER=2
while [ $ITER -lt 10 ]; do
	$BWA index assembly.polished$ITER.fasta
	$BWA mem -t $CPUS assembly.polished$ITER.fasta F.cor.fq R.cor.fq | $SAMTOOLS sort -o alignment.sort.bam
	$SAMTOOLS index alignment.sort.bam
	$PILON --genome assembly.polished$ITER.fasta --frags alignment.sort.bam --output assembly --outdir pilonOut --changes --fix snps,indels,gaps --iupac --threads $CPUS
	cat pilonOut/assembly.changes | nodejs ../pilonChanges.js >> polish.changes.txt
	ITER=`expr $ITER + 1`
	sed -i 's/_pilon//g' pilonOut/assembly.fasta
	cp pilonOut/assembly.fasta assembly.polished$ITER.fasta
done
cp assembly.polished$ITER.fasta assembly.polished.fasta
rm alignment.sort.bam
rm alignment.sort.bam.bai
rm assembly.polish*.fasta.amb
rm assembly.polish*.fasta.ann
rm assembly.polish*.fasta.bwt
rm assembly.polish*.fasta.pac
rm assembly.polish*.fasta.sa

#Final polish for non-IUPAC output and hapolid variants
$BWA index assembly.polished.fasta
$BWA mem -t $CPUS assembly.polished.fasta F.cor.fq R.cor.fq | $SAMTOOLS sort -o alignment.sort.bam
$SAMTOOLS index alignment.sort.bam
$PILON --genome assembly.polished.fasta --frags alignment.sort.bam --output assembly --outdir pilonOut --changes --fix snps,indels,gaps --threads $CPUS
echo "Non-IUPAC final polish" >> polish.changes.txt
cat pilonOut/assembly.changes | nodejs ../pilonChanges.js >> polish.changes.txt
cp pilonOut/assembly.fasta assembly.final.fasta
$PILON --genome assembly.polished.fasta --frags alignment.sort.bam --output assembly --outdir pilonOut --variant --vcf --diploid --threads $CPUS
$SAMTOOLS depth -aa alignment.sort.bam | nodejs coverage.js > coverage.txt
