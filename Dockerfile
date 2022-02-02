FROM python:3.8.6-buster

MAINTAINER ken-ichi chiba <kchiba@hgc.jp>

RUN apt-get -y update
RUN apt-get install -y liblzma-dev libbz2-dev ncurses-dev zlib1g-dev
RUN apt-get install -y wget 

RUN mkdir -p /tools

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install pysam==0.15.4

RUN cd /tools && \
echo "build samtools-1.9" && \
wget -nc https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2 && \
tar xjvf samtools-1.9.tar.bz2 && \
rm -f samtools-1.9.tar.bz2 && \
cd samtools-1.9 && \
make && \
make install

RUN cd /tools && \
echo "build htslib-1.9" && \
wget -nc https://sourceforge.net/projects/samtools/files/samtools/1.9/htslib-1.9.tar.bz2 && \
tar xvf htslib-1.9.tar.bz2 && \
rm -f htslib-1.9.tar.bz2 && \
cd htslib-1.9 && \
make && \
make install

RUN cd /tools && \
echo "build bcftools-1.9" && \
wget -nc https://sourceforge.net/projects/samtools/files/samtools/1.9/bcftools-1.9.tar.bz2 && \
tar xvf bcftools-1.9.tar.bz2 && \
rm -f bcftools-1.9.tar.bz2 && \
cd bcftools-1.9 && \
./configure && \
make && \
make install

RUN cd /tools && \
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
tar -zxvf bedtools-2.29.1.tar.gz && \
cd bedtools2 && \
make
ENV PATH $PATH:/tools/bedtools2/bin

RUN cd /tools && \
echo "build ob_utils-0.0.11" && \
wget -nc https://github.com/ken0-1n/ob_utils/archive/v0.0.11.tar.gz && \
tar xzvf v0.0.11.tar.gz && \
cd ob_utils-0.0.11 && \
python3 setup.py install

