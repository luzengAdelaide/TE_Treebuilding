# Example of building tree by using CARP, job running on HPC

# Build L2 tree based on genomic 
# Use echinda genome as an example

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=10GB

# Notification configuration
#SBATCH --mail-type=END                                        
#SBATCH --mail-type=FAIL                                       
#SBATCH --mail-user=lu.zeng@adelaide.edu.au  

# Loadsoftwares need to use later
module load MUSCLE
module load HMMER
module load BEDTools
module load CENSOR
module load wu-blast
module load bioperl


# Extract L2 sequences that is longer 2kb from the Echidna repeat library
# Extract all L2 sequences from CARP output
perl extract_l2_repeat_fasta.pl Echidna_len85_Library.fasta_repeat > L2_echidna.fa
# Extract L2 sequences that are longer than 2kb
perl extract_2kb_sequence.pl L2_echidna.fa
# Rename the long L2 sequences for next step
perl extract_family_info.sh 2kb.L2_echidna.fa > title_L2_echidna
# Extract the original genome sequences from L2 consensus sequences
rsync --files-from=title_L2_echidna /genomes/echidna/mfa/ ./
# Merge them together
cat *.mfa > 2kb.L2_echidna_genome.fa
# Reverse minus strand
# Use CENSOR to figure out the strand direction, use crocodile consensus sequences as index
censor -bprm cpus=8 -lib L2_echidna.fa 2kb.L2_echidna_genome.fa
perl divide-+.pl 2kb.L2_echidna_genome.fa.map
perl extract_fasta_seq.pl 2kb.L2_echidna_genome.fa.map_comple 2kb.L2_echidna_genome.fa > comple_tachy_l2 
perl extract_fasta_seq.pl 2kb.L2_echidna_genome.fa.map_forward 2kb.L2_echidna_genome.fa > forward_tachy_l2
perl ../reverse_complement.pl comple_tachy_l2 > reverse_tachy_l2

# Add outgroup from Repbase L2 (e.g. Crocodile)
# cat repbase_L2_tree reverse_tachy_l2 forward_tachy_l2 > 2kb.L2_echidna_genome.use

# Merge echidna L2 sequences together
cat reverse_tachy_l2 forward_tachy_l2 > 2kb.L2_echidna_genome.use
# Use muscle to run alignments
muscle -in 2kb.L2_echidna_genome.use -out 2kb.L2_echidna.afa -maxiters 2

# Trim the alignment outputs
sed -i 's/\s.*$//' 2kb.L2_echidna.afa
Gblocks 2kb.L2_echidna.afa -t=d -p=n -e=.gb -b5=h
tr -d " \t" < 2kb.L2_echidna.afa.gb > 2kb.L2_echidna.gbHalf.afa

# Build phylogeny tree
fasttree -gtr -nt 2kb.L2_echidna.gbHalf.afa > gblock_tree_L2_denovo

# Extract ORF2 sequences
usearch -fastx_findorfs 2kb.L2_echidna_genome.use -aaout tachy.ORF2.candidate.translated -mincodons 500
hmmscan --domtblout tachy.ORF2.candidate.translated.out /data/rc003/lu/pfam/Pfam-A.hmm tachy.ORF2.candidate.translated > tachy.ORF2.candidate.translated.log
sed -i '/^#/ d' tachy.ORF2.candidate.translated.out
awk '{print $1 "\t" $4}' tachy.ORF2.candidate.translated.out > tachy_ORF2_filtered.tmp
awk -v i=2 'NR>1 && $i!=p { print "" }{ p = $i } 1' tachy_ORF2_filtered.tmp > tachy_ORF2_sep.tmp
cat tachy_ORF2_sep.tmp | cut -f 1,2| awk '{if ( $1=="RVT_1" || $1=="RVT_3" ) print $2 }' | sort -u > tachy_confirmedORF2.txt
cat tachy.ORF2.candidate.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "tachy_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > tachy_confirmedORF2.fasta

# Extract RT Domain
hmmscan --domtblout tachy_confirmedORF2.fasta.out /data/rc003/lu/pfam/Pfam-A.hmm tachy_confirmedORF2.fasta > tachy_confirmedORF2.fasta.log
cat tachy_confirmedORF2.fasta.out | awk '{if (($1=="RVT_1" || $1=="RVT_3") && ($21-$20>200)) print $4 "\t"  $20 "\t" $21}' > tachy_RT_domain.bed
fastaFromBed -fi tachy_confirmedORF2.fasta -bed tachy_RT_domain.bed -fo tachy_RT_domain.fasta

# Extract the sequence name from RT domain output
grep '>' tachy_RT_domain.fasta > title_RTdomain
sed -i 's/|.*$//; s/>//' title_RTdomain 
tr '\n' ',' < title_RTdomain > title_RTdomain.use
