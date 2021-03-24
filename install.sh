#!/bin/bash


projdir=`pwd`

bindir=$projdir/src/caus-bin/
mkdir -p $bindir
mkdir -p $projdir/src/log/

errs=0

##### Download and install minimap2 ######

echo "Downloading and installing minimap2"
if [[ ! -s $bindir/minimap2 ]]; then

    if [[ ! -d $projdir/src/minimap2 ]]; then
	cd $projdir/src/
	git clone https://github.com/lh3/minimap2.git &> $projdir/src/log/minimap2_cloning.log
    fi

    if [[ ! -s $projdir/src/minimap2/minimap2 ]]; then
	cd $projdir/src/minimap2
	make &> $projdir/src/log/minimap2_installation.log
    fi

    chmod 755 minimap2
    cp minimap2 $bindir
fi

if  [[ ! -s $bindir/minimap2 ]]; then
    echo " !! Error: minimap2 not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if minimap2 was downloaded properly:" $projdir/src/log/minimap2_cloning.log 
    echo "   Check if the minimap2 was compiled properly:" $projdir/src/log/minimap2_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/minimap2/minimap2 $bindir/minimap2 
    
    errs=$(($errs+1))
else
    echo " minimap2 succesfully installed!"
    rm -rf $projdir/src/minimap2/
fi

##### Download and install SMALT ######
echo; echo "Downloading and installing Smalt"
if [[ ! -s $bindir/smalt ]]; then
   
    if [[ ! -d $projdir/src/smalt-0.7.4 ]]; then
        cd $projdir/src/
        wget ftp://ftp.sanger.ac.uk/pub/resources/software/smalt/smalt-0.7.4.tgz &> $projdir/src/log/smalt_wget.log
        tar -xvzf smalt-0.7.4.tgz &> $projdir/src/log/smalt_untar.log
        rm -f smalt-0.7.4.tgz
    fi

    cp $projdir/src/smalt-0.7.4/smalt_x86_64 $bindir/smalt
    chmod 755 $bindir/smalt
fi
if  [[ ! -s $bindir/smalt ]]; then 
    echo " !! Error: smalt not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if smalt was downloaded properly:" $projdir/src/log/smalt_wget.log 
    echo "   Check if the folder was uncompressed properly:" $projdir/src/log/smalt_untar.log

    # Cleaning up       
    rm -rf $projdir/src/smalt-0.7.4/ $bindir/smalt

    errs=$(($errs+1))
else
    echo " Smalt succesfully installed!"
    rm -rf $projdir/src/smalt-0.7.4/
fi


###### Compile caus sources ######

echo; echo "Compiling caus sources"

srcs=( caus_seqout fastq2fasta caus_fasta caus_shred caus_assign caus_smalt0 caus_smalt1 caus_smalt2 caus_smalt3 caus_clean caus )

cd $projdir/src
make &> $projdir/src/log/sources_compilation.log

echo; echo "Checking installation:"
for src in "${srcs[@]}"; do
    if [[ ! -s $bindir/$src ]]; then 
        echo " !! Error: executable $src missing in $bindir"
	echo "    Please check for errors the log file:" $projdir/src/log/sources_*	
        errs=$(($errs+1))
    fi
done

if [  $errs -gt 0 ]; then echo; echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi




