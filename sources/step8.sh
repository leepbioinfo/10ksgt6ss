#!/bin/bash

prefix=20240211 
wd=$(readlink -f $(pwd))
bd=$(basename $wd)
lib=$bd/data/models

# Setup libraries at $bd/data/models/
#
# The klist of model names is available as text files in ../data/models
# Pfam, CDD and TXSScan (rocha) models are available from their repositories
# Aravind models were kindly provided by Dr. L. Aravind Iyer from NCBI
#
for list in $lib/*.lst
do
	library=$(basename $list .lst)
	cd $lib
	mkdir $library
	cd $library
	if [ "$library" == "st" ]; then
		ln -s ../../../hmm/*.hmm .
	fi
	cd $lib
	cat $library/*.hmm > ${library}.hmm
	hmmpress ${library}.hmm
done

# Run all models
cd $wd
for library in pfam aravind rocha cdd st
do
	mkdir $library
	cd $library
	\ls $lib/$library/*.hmm \
	| parallel -N1 -j36 '
	   bn={/.}.merged2.hmmsearch;
	   if [ -f ${bn}.out ]; then
		echo $bn.out done;
	   else
		echo Processing $bn;
		hmmsearch -o ${bn}.out {} ../../data/merged2.c100i100.fa;
		hmmer2table -y ${bn}.out > ${bn}.tsv;
	   fi
	'
	find . -name '*.hmmsearch.tsv' -size 0c -exec rm -f \{\} \;
	cd $wd
done

# Collect HMMsearch results
cd $wd
awk '!/^sequence/ && $5<=0.01{print $1}' {aravind,cdd,pfam,rocha,st}/*.hmmsearch.tsv | sort -u > ${prefix}.acc
if [ ! -f $bd/data/merged2.faa.ssi ]; then
  esl-sfetch --index $bd/data/merged2.faa
fi
cat ${prefix}.acc \
| parallel -N1 -j36 esl-sfetch $bd/data/merged2.faa \
> ${prefix}.combined.fa 2> ${prefix}.combined.eslsfetch.err
pfetch $(awk '/^seq/{print $2}' ${prefix}.combined.eslsfetch.err) >> ${prefix}.combined.fa

# Run scanseqs to identify best matching models
mkdir scanseqs
cd scanseqs
scanseqs -s none -s phobius=phobius \
 -s aravind=hmmscan:$lib/aravindDB.hmm \
 -s aravind=rpsblast:$lib/aravindDB \
 -s pfam=hmmscan:$lib/Pfam-A.hmm \
 -s cdd=hmmscan:$lib/cdd.hmm \
 -s rocha=hmmscan:$lib/rocha.hmm \
 -s st=hmmscan:$lib/st.hmm \
 $prefix ${prefix}.combined.fa

exit $?
