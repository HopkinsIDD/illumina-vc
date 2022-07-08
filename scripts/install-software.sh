#!/bin/bash

## For installing software needed for genome assembly from Illumina reads
## Installed software will only be checked for updating


# ---------------- Input Paramters and Files  ---------------- #

DIR=$1 # directory into which software should be installed



# check that an installation directory is provided
if [ $# -eq 0 ]; then
  echo "ERROR: no install directory supplied"
  exit
fi


# remove the trailing slash from the directory name if there is one
DIR=${DIR%/}


# make the installation directory if it does not yet exist
if [ ! -d $DIR ]; then
  mkdir -p $DIR
fi

### DEFINE FUNCTIONS ###

# function to install linux packages not already installed
install_if_not_installed(){
  echo $1
  apt -qq list $1 2> /dev/null | (grep -qE "installed" ) || apt install $1
}

# function to initiate software installation from github
enter_or_clone(){
  echo $1/$2
  cd $DIR
  if [[ ! -d $2 ]] ; then
    git clone --recursive https://github.com/$1/$2.git
    cd $2
  else
    cd $2
  fi
}

# function to fully install software from github
pull_and_install(){
  echo $1
  cd $DIR/$1
  rm -f .initialpull
  if [[ ! $(git pull | grep "Already up to date") ]] || [[ ! -e .initialpull ]] ; then 
    if [[ ! -e .initialpull ]]; then
      for cmd in "${@:2}"; do
        $cmd
      done
    fi
  else
    echo "  Up to date"
  fi
  touch .initialpull
}


### INSTALL BASIC LINUX AND PYTHON PACKAGES ###

echo "installing linux packages with apt"

# get a list of already-installed packages (Necessary?)
cat /etc/apt/sources.list | sed "s/restricted/restricted universe/" | sed "s/universe universe/universe/" &> /etc/apt/sources.list.new
mv /etc/apt/sources.list.new /etc/apt/sources.list

# list of packages to install
PACKAGES=(bioperl cmake curl default-jre emboss gcc git gzip g++ libbz2-dev libcurl4-gnutls-dev libffi-dev libfile-slurp-perl libjson-perl liblwp-protocol-https-perl liblzma-dev libncurses5-dev libssl-dev libtext-csv-perl libwww-perl libz-dev libeigen3-dev make ncbi-blast+ python-virtualenv python3-all-dev python3-pip unzip wget zlib1g snp-sites mauve figtree autoconf)

# update apt and install packages with apt
apt update
for PKG in ${PACKAGES[@]}; do
  install_if_not_installed $PKG
done

echo "installing python packages"

# intall python packages with pip
sudo -H -u Illumina python3 -m pip install --upgrade --user pip
sudo -H -u Illumina python3 -m pip install --upgrade --user setuptools
PYTHON_PACKAGES=(numpy pandas biopython pysam)

for PKG in ${PYTHON_PACKAGES[@]}; do
  sudo -H -u Illumina python3 -m pip install --user $PKG
done


### INSTALL SOFTWARE FROM GITHUB ###

## # Install latest bwa
enter_or_clone lh3 bwa
pull_and_install bwa make

## # Install latest htslib
enter_or_clone samtools htslib
pull_and_install htslib autoheader autoconf ./configure make "make install"

## # Install latest samtools
enter_or_clone samtools samtools
pull_and_install samtools autoheader "autoconf -Wno-syntax" ./configure make "make install"

## # Install latest picard  
enter_or_clone broadinstitute picard
pull_and_install picard "gradlew shadowJar"

## # Install bcftools
enter_or_clone samtools bcftools
pull_and_install bcftools "autoreconf -i -Wno-syntax" ./configure make "make install"

## # Install Snippy (Should we do this to all softwares? Or should we add the software directory in search path)
enter_or_clone tseemann snippy
pull_and_install snippy "cp -R * /usr/local/"

## # Install Bedtools
cd $DIR
wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
tar -xzvf bedtools-2.28.0.tar.gz
cd bedtools2
make

## # Install RAxML (? Version)
enter_or_clone stamatak standard-RAxML
pull_and_install standard-RAxML "make -f Makefile.gcc" "rm *.o"


## # Install gubbins
enter_or_clone sanger-pathogens gubbins
pull_and_install gubbins "autoreconf -i" "./configure" "make" "make install" "cd python" "python3 setup.py install"

## # Install IQ-TREE
enter_or_clone Cibiv IQ-TREE
pull_and_install IQ-TREE "mkdir build" "cd build" "cmake -DIQTREE_FLAGS=omp .." "make -j4" "make install"

## # Install FigTree
enter_or_clone rambaut figtree
pull_and_install figtree "cd release"

echo "installation complete; check above output for possible errors"