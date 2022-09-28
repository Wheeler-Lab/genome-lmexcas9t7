set -e

mkdir work || echo "Work directory exists"
cd work

python3 -m venv venv

#For nanopore qc
./venv/bin/python3 -m pip install git+https://github.com/rrwick/Porechop.git

#For minimap/miniasm build
git clone https://github.com/lh3/minimap2
cd minimap2 && make && cd ..
git clone https://github.com/lh3/miniasm
cd miniasm && make && cd ..

#For racon consensus
git clone --recursive https://github.com/lbcb-sci/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ../..

#For pilon polishing
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

./venv/bin/python3 -m pip install cutadapt

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


# Versions
minimap2/minimap2 -V > version.minimap2.txt || echo
miniasm/miniasm -V > version.miniasm.txt || echo
racon/build/bin/racon --version > version.racon.txt || echo
./bwa/bwa 2>&1 ./bwa/bwa 2>&1 | grep Version > version.bwa.txt || echo
./samtools/bin/samtools 2>&1 | grep Version > version.samtools.txt || echo
java -jar pilon/pilon.jar | grep version > version.pilon.txt || echo
./venv/bin/porechop --version > version.porechop.txt || echo
