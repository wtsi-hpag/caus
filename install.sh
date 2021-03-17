#!/bin/bash


projdir=`pwd`

bindir=$projdir/caus-bin/
mkdir -p $bindir
mkdir -p $projdir/log/

errs=0

##### Download and install minimap2 ######

echo "Downloading and installing minimap2 "
if [[ ! -s $bindir/minimap2 ]]; then

    if [[ ! -d $projdir/minimap2 ]]; then
	cd $projdir/
	git clone https://github.com/lh3/minimap2.git &> $projdir/log/minimap2_cloning.log
    fi

    if [[ ! -s $projdir/minimap2/minimap2 ]]; then
	cd $projdir/minimap2
	make &> $projdir/log/minimap2_installation.log
    fi

    cp minimap2 $bindir
fi

if  [[ ! -s $bindir/minimap2 ]]; then
    echo " !! Error: minimap2 not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if minimap2 was downloaded properly:" $projdir/log/minimap2_cloning.log 
    echo "   Check if the minimap2 was compiled properly:" $projdir/log/minimap2_installation.log

    # Cleaning up
    cd $projdir/
    rm -rf $projdir/minimap2/minimap2 $bindir/minimap2 
    
    errs=$(($errs+1))
else
    echo " minimap2 succesfully installed!"
    rm -rf $projdir/minimap2/
fi


##### Download and install SMALT ######
echo; echo "Downloading and installing Smalt"
if [[ ! -s $bindir/smalt ]]; then
   
    if [[ ! -d $projdir/smalt-0.7.4 ]]; then
	cd $projdir/
	wget ftp://ftp.sanger.ac.uk/pub/resources/software/smalt/smalt-0.7.4.tgz &> $projdir/log/smalt_wget.log
	tar -xvzf smalt-0.7.4.tgz &> $projdir/log/smalt_untar.log
	rm -f smalt-0.7.4.tgz
    fi

    cp $projdir/smalt-0.7.4/smalt_x86_64 $bindir/smalt
fi
if  [[ ! -s $bindir/smalt ]]; then 
    echo " !! Error: smalt not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if smalt was downloaded properly:" $projdir/log/smalt_wget.log 
    echo "   Check if the folder was uncompressed properly:" $projdir/log/smalt_untar.log

    # Cleaning up	
    rm -rf $projdir/smalt-0.7.4/ $bindir/smalt

    errs=$(($errs+1))
else
    echo " Smalt succesfully installed!"
    rm -rf $projdir/smalt-0.7.4/
fi


###### Compile CAUS sources ######

echo; echo "Compiling caus sources"

srcs=(caus_seqout caus_fasta caus_shred caus_assign caus_smalt0 caus_smalt1 caus_smalt2 caus_smalt3 caus_clean caus)

cd $projdir/
make &> $projdir/log/sources_compilation.log





