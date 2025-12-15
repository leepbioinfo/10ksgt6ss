#!/bin/bash
# File: /projects/salmonella/work/20210212/source1.sh

# Basename
target=merged2.c100i100
input=${target}.fa

# Detecting sequence signatures in the new proteins

# Phobius
cat ${input} | parallel -N1 -j10 --pipe --recstart ">" phobius > ${target}.phobius.out 2> ${target}.phobius.err
phobius2table -e 0.0101 ${target}.phobius.out > ${target}.phobius.tsv

# Pfam: hmmscan
cat ${target}.fa \
| parallel -N1 -j10 --pipe --recstart ">" hmmscan /databases/pfam/Pfam-A.hmm - \
> ${target}.pfam.hmmscan.out \
2> ${target}.pfam.hmmscan.err
hmmer2table ${target}.pfam.hmmscan.out > ${target}.pfam.hmmscan.tsv
domain2architecture -e 0.0101 ${target}.phobius.tsv ${target}.pfam.hmmscan.tsv > ${target}.pfam.hmmscan.arch
architecture2table ${target}.pfam.hmmscan.arch > ${target}.pfam.hmmscan.arch.tsv
ln -s ${target}.pfam.hmmscan.arch ${target}.pfam.scan.arch
ln -s ${target}.pfam.hmmscan.arch.tsv ${target}.pfam.scan.arch.tsv

# Aravind: hmmscan
cat ${target}.fa \
| parallel --pipe -N1 -j36 --recstart '>' hmmscan --cpu 1 /databases/profiledb/hmmer/aravinddb - \
> ${target}.aravind.hmmscan.out \
2> ${target}.aravind.hmmscan.err
hmmer2table ${target}.aravind.hmmscan.out > ${target}.aravind.hmmscan.tsv

# Aravind: rpsblast
cat ${target}.fa \
| parallel -N1 -j36 --pipe --recstart '>' rpsblast -db /databases/profiledb/rpsdb/aravinddb \
> ${target}.aravind.rpsblast.out \
2> ${target}.aravind.rpsblast.err
cat ${target}.aravind.rpsblast.out | parallel -N1 -j36 --pipe --recstart RPSBLAST blast2table -s > ${target}.aravind.rpsblast.tsv

# Aravind: domain2architecture
(cat ${target}.phobius.tsv; awk '$5<=0.1{print}' ${target}.aravind.rpsblast.tsv) \
| domain2architecture -e 0.0101 > ${target}.aravind.scan.arch
architecture2table ${target}.aravind.scan.arch > ${target}.aravind.scan.arch.tsv

# Search for matches to Rocha/TXSScan T6SS models from:
# https://github.com/macsy-models/TXSScan
\ls -1 /databases/macsydata/TXSScan/profiles/T6SSi*.hmm | parallel -N1 -j20 "
 hmmsearch -o ${target}.rocha.{/.}.hmmsearch.out {} ${input} 2> ${target}.rocha.{/.}.hmmsearch.err;
 hmmer2table -r {/.} ${target}.rocha.{/.}.hmmsearch.out > ${target}.rocha.{/.}.hmmsearch.tsv
"
domain2architecture -e 0.1 ${target}.rocha.*.hmmsearch.tsv > ${target}.rocha.hmmsearch.arch
architecture2table ${target}.rocha.hmmsearch.arch > ${target}.rocha.hmmsearch.arch.tsv
ln -s ${target}.rocha.hmmsearch.arch ${target}.rocha.scan.arch
ln -s ${target}.rocha.hmmsearch.arch.tsv ${target}.rocha.scan.arch.tsv

# Run remaining Pfam models
cat ../../data/models/pfam.lst | parallel -N1 -j20 "
 hmmfetch /databases/pfam/Pfam {} | hmmsearch --cpu 1 -o ${target}.pfam.{/.}.hmmsearch.out - ${input} 2> ${target}.pfam.{/.}.hmmsearch.err;
"
hmmer2table ${target}.pfam.*.hmmsearch.out > ${target}.pfam.hmmsearch.tsv
domain2architecture -e 0.1 ${target}.pfam.hmmsearch.tsv > ${target}.pfam.hmmsearch.arch

# Searching homologs of CDD MIX and FIX models
tfilter -f '$F[2]=~/MIX|FIX-like/' -f 'use File::Basename;$F[1]=dirname($H{"current_file"})."/$F[1].hmm"' ../data/models/cdd/cddid_all.t6ss.tsv -i 1 \
| parallel -N1 hmmsearch -o ${target}.cdd.{/.}.hmmsearch.out {} ${input}
hmmer2table ${target}.cdd.*.hmmsearch.out > ${target}.cdd.hmmsearch.tsv
domain2architecture -e 1 ${target}.cdd.hmmsearch.tsv > ${target}.cdd.hmmsearch.arch
architecture2table ${target}.cdd.hmmsearch.arch > ${target}.cdd.hmmsearch.arch.tsv
ln -s ${target}.cdd.hmmsearch.arch ${target}.cdd.scan.arch
ln -s ${target}.cdd.hmmsearch.arch.tsv ${target}.cdd.scan.arch.tsv

# Aravind's models for PAAR and VgrG are the same as the ones in Pfam
# I could not find any other T6SS signatures

# Merge results from all searches (see source1.sh for Rocha and source2.sh for Pfam)
cut -f 1 ${target}.*.hmmsearch.arch | sort -u | grep -v ^ID > ../data/t6ss.acc

# Remove empty files
find . -size 0c -exec rm -f \{\} \;
