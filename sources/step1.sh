#!/bin/bash
# File: /projects/salmonella/data/merged2/source1.sh

# Extracting data from GFF files of the initial set of genomes
# 1. 10.000 Salmonella Genomes Project (10KSG), provided by Jay Hilton
# 2. S. Enteritidis genomes assembled and annotated by Gemma C. Langridge et. al. (PMID: 25535353)
# 3. RefSeq Genomes GCF_000006945.2, GCF_000022165.1 and GCF_000210855.2
#
# !!!!! IMPORTANT NOTE !!!!!
# Datsets 2 and 3 where NOT included in our genome analysis or used to 
# identify new toxin domains but where not removed from merged2.tsv and
# ssg.tsv. Therefore, proteins from these genomes may appear as cluster
# representatives but such proteins where replaced by identical
# proteins from the 10KSG projectt

# Parse
\ls -1 ../data/*.gff | parallel -N1 -j30 rgaparser --assembly $(echo {/.} | sed "s/_genomic//") --features {/.}.a2o.tsv --dna {/.}.fna --protein {/.}.faa {}
\ls -1 ../data/*.gbff | parallel -N1 -j30 rgaparser --informat genbank --assembly $(echo {/.} | sed "s/_genomic//") --features {/.}.a2o.tsv --dna {/.}.fna --protein {/.}.faa {}

# Concatenate features tables
(head -n1 GCF_000210855.2_ASM21085v2_genomic.a2o.tsv; cat *.a2o.tsv | grep -v ^nucleotide) > merged2.tsv

# Merging protein FASTA files
cat *.faa > merged2.faa

# Cluster sequences using hashclust and 100% coverage + 100% identity
# Table mapping all sequences to non-redundant representatives is mergednrdb.tsv
# Non-redundant FASTA file is mergednrdb_rep.fasta
mmseqs createdb merged2.faa merged2
mmseqs clusthash merged2 merged2nr --threads 20 --min-seq-id 1
mmseqs clust merged2 merged2nr merged2.c100i100 --threads 20 --cluster-mode 1 --max-iterations 1000000
mmseqs createtsv merged2 merged2 merged2.c100i100 merged2.c100i100.tsv
mmseqs createsubdb merged2.c100i100 merged2 merged2.c100i100.seqs
mmseqs convert2fasta merged2.c100i100.seqs merged2.c100i100.fa
mmseqs easy-cluster merged2.c100i100.fa merged2.c80e3 tmp --threads 36 > merged2.c80e3.log 2>&1 
ln -s merged2.c80e3_rep_seq.fasta merged2.c80e3.fa
ln -s merged2.c80e3_cluster.tsv merged2.c80e3.tsv
rm -fr tmp
