set -e

mkdir work
cd work

python3 -m venv ./venv
./venv/bin/python3 -m pip install pandas xlrd openpyxl

#Fetch TriTrypDB v6.0
#Old genome version (LmxM.00 contig)
mkdir RefGFF_6.0
curl https://tritrypdb.org/common/downloads/release-6.0/LmexicanaMHOMGT2001U1103/fasta/data/TriTrypDB-6.0_LmexicanaMHOMGT2001U1103_Genome.fasta -o RefGFF_6.0/genome.fasta
curl https://tritrypdb.org/common/downloads/release-6.0/LmexicanaMHOMGT2001U1103/gff/data/TriTrypDB-6.0_LmexicanaMHOMGT2001U1103.gff -o RefGFF_6.0/genome.gff

#Fetch TriTrypDB v50
#New genome version (mexFOS1_* contigs)
mkdir RefGFF_50
curl https://tritrypdb.org/common/downloads/release-50/LmexicanaMHOMGT2001U1103/fasta/data/TriTrypDB-50_LmexicanaMHOMGT2001U1103_Genome.fasta -o RefGFF_50/genome.fasta
curl https://tritrypdb.org/common/downloads/release-50/LmexicanaMHOMGT2001U1103/gff/data/TriTrypDB-50_LmexicanaMHOMGT2001U1103.gff -o RefGFF_50/genome.gff

#Determine coordinate remapping for two reference databases
nodejs ../determineGFFRemap.js RefGFF_6.0/genome.gff RefGFF_50/genome.gff > coordRemap.json

#Fetch Fiebig et al. data
mkdir FiebigGFF
curl https://doi.org/10.1371/journal.ppat.1005186.s005 -L -o FiebigGFF/s5.slas.xlsx
curl https://doi.org/10.1371/journal.ppat.1005186.s006 -L -o FiebigGFF/s6.pas.xlsx
curl https://doi.org/10.1371/journal.ppat.1005186.s007 -L -o FiebigGFF/s7.extended.xlsx
curl https://doi.org/10.1371/journal.ppat.1005186.s008 -L -o FiebigGFF/s8.truncated.xlsx
curl https://doi.org/10.1371/journal.ppat.1005186.s009 -L -o FiebigGFF/s9.novel.xlsx
cd FiebigGFF
for f in *.xlsx
do
    ../venv/bin/python3 -c "import pandas; df = pandas.read_excel('$f'); df.to_csv('$f.csv', index=False);"
done
cd ..

#Correct/tidy Fiebig et al. data to correctly formatted GFFs
nodejs ../correctSlasPas.js FiebigGFF/s5.slas.xlsx.csv > FiebigGFF/s5.slas.gff
nodejs ../correctSlasPas.js FiebigGFF/s6.pas.xlsx.csv > FiebigGFF/s6.pas.gff
nodejs ../correctCds.js FiebigGFF/s7.extended.xlsx.csv > FiebigGFF/s7.extended.gff
nodejs ../correctCds.js FiebigGFF/s8.truncated.xlsx.csv > FiebigGFF/s8.truncated.gff
nodejs ../correctNovel.js FiebigGFF/s9.novel.xlsx.csv > FiebigGFF/s9.novel.gff

#Remap Fiebig GFFs to new reference coordinates
cd FiebigGFF
for f in *.gff
do
	nodejs ../../remapGFFCoords.js $f ../coordRemap.json > $f.cgff
done
cd ..

#Parse GFFs to a json format
nodejs ../parseGFF.js RefGFF_6.0/genome.gff > RefGFF_6.0/genome.gff.json
nodejs ../parseGFF.js RefGFF_50/genome.gff > RefGFF_50/genome.gff.json
cd FiebigGFF
for f in *.cgff
do
	nodejs ../../parseGFF.js $f > $f.json
done
cd ..

#Merge novel CDS, PAS and SLAS entries with reference
nodejs ../mergeJSON.js RefGFF_50/genome.gff.json FiebigGFF/s9.novel.gff.cgff.json > genome.1.tmp.json
nodejs ../mergeJSON.js genome.1.tmp.json FiebigGFF/s5.slas.gff.cgff.json > genome.2.tmp.json
nodejs ../mergeJSON.js genome.2.tmp.json FiebigGFF/s6.pas.gff.cgff.json > genome.3.tmp.json

#Modify CDSs, if new CDS is an orf
nodejs ../updateCDS.js genome.3.tmp.json FiebigGFF/s7.extended.gff.cgff.json RefGFF_50/genome.fasta > genome.4.tmp.json
nodejs ../updateCDS.js genome.4.tmp.json FiebigGFF/s8.truncated.gff.cgff.json RefGFF_50/genome.fasta > genome.5.tmp.json

#Correct gene coordinates and add UTRs
nodejs ../updateCoords.js genome.5.tmp.json > genome.6.tmp.json

#Convert json to GFF
nodejs ../outputGFF genome.6.tmp.json > genome.gff
grep "gene" genome.gff > genome.gene.gff
grep "CDS" genome.gff > genome.cds.gff
grep "PAS" genome.gff > genome.pas.gff
grep "SLAS" genome.gff > genome.slas.gff
grep "_prime_UTR" genome.gff > genome.utr.gff

#Tidy JSONs
mv genome.6.tmp.json genome.json
rm *.tmp.json
