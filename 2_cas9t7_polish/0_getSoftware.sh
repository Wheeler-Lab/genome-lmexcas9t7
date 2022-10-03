set -e

mkdir work
cd work

mkdir pilon
cd pilon && curl https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar -L -o pilon.jar && cd ..

curl https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz -L -o salmon.tar.gz
tar -xf salmon.tar.gz

git clone https://github.com/mourisl/rcorrector.git
cd rcorrector
make
cd ..

curl https://github.com/FelixKrueger/TrimGalore/archive/0.6.0.zip -L -o trimgalore.zip
unzip trimgalore.zip

python3 -m venv ./cutadapt_venv
./cutadapt_venv/bin/python3 -m pip install cutadapt

git clone https://github.com/lh3/bwa.git
cd bwa; make
cd ..

curl https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 -L -o samtools.tar.bz2
tar -xf samtools.tar.bz2
cd samtools-1.16
./configure --prefix=`pwd`/../samtools
make
make install
cd ..

git clone https://github.com/harvardinformatics/TranscriptomeAssemblyTools

git clone https://github.com/ncbi/ncbi-vdb
git clone https://github.com/ncbi/sra-tools
cd ncbi-vdb
./configure --relative-build-out-dir
make -j12
cd ../sra-tools
./configure --relative-build-out-dir
make -j12
make
cd ..
./OUTDIR/sra-tools/linux/gcc/x86_64/rel/bin/vdb-config > /dev/null  # Trick SRA-TOOLS into thinking we've configured it...