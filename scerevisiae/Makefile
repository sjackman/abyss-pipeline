# Number of threads
j=12

# ABySS parameters
k=512
K=96
e=2
E=0
c=1.01
l=40
s=500
n=10
S=500
N=100

# S288C reference genome
ref=R64

# Targets

all: \
	assembly-stats.tsv \
	assembly-stats.html \
	scerevisiae.quast/transposed_report.tsv \
	quast-stats.html

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Download data

# Download the reference genome from Ensembl
$(ref).fa:
	curl ftp://ftp.ensembl.org/pub/release-78/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz |gunzip -c >$@

# Download gene annotations from Ensembl in GTF format
$(ref).gtf:
	curl ftp://ftp.ensembl.org/pub/release-78/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.78.gtf.gz |gunzip -c >$@

# Download gene annotations from Saccharomyces Genome Database (SGD)
$(ref).gff:
	curl -o $@ http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff

# Download Nanocorr corrected reads
data/W303_ONT_Nanocorr_Corrected_reads.fa.gz:
	curl -o $@ http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/W303_ONT_Nanocorr_Corrected_reads.fa.gz

# Download Nanocorr CA spec file
w303_ont.spec:
	curl -o $@ http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/w303_ont.spec

# Download Nanocorr polished assembly
data/W303_ONT_Assembly_CA_polished.fa.gz:
	curl -o $@ http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/W303_ONT_Assembly_CA_polished.fa.gz

# BWA

# Index the reference
%.fa.bwt: %.fa
	bwa index $<

# Align reads to the reference using BWA MEM
$(ref).%.sam: data/%.fa.gz $(ref).fa.bwt
	bwa mem -t$j $(ref).fa $< >$@

# Sort a SAM file and create a BAM file
%.sort.bam: %.sam
	samtools view -Su $< |samtools sort - $*.sort

# Index a BAM file
%.bam.bai: %.bam
	samtools index $<

# bedtools

# Determine regions with no coverage
%.bam.depth0.bed: %.bam
	bedtools genomecov -bga -ibam $< |awk '$$4 == 0' |cut -f1-3 >$@

# Determine regions with 1+ coverage
%.bam.depth1+.bed: %.bam
	bedtools genomecov -bga -ibam $< |awk '$$4 > 0' |cut -f1-3 |bedtools merge >$@

# Determine regions with 2+ coverage
%.bam.depth2+.bed: %.bam
	bedtools genomecov -bga -ibam $< |awk '$$4 > 1' |cut -f1-3 |bedtools merge >$@

# Extract FASTA sequences from BED regions
%.bed.fa: %.bed
	bedtools getfasta -fi $(ref).fa -bed $< -fo $@

# Determine the mean BAM coverage of BED regions
%.bam.coverage: %.bam %.bam.depth1+.bed
	bedtools coverage -d -abam $< -b $*.bam.depth1+.bed |datamash --header-out -g 1,2,3 mean 5 >$@

# Correct reads using BFC (Bloom filter corrector)
%.bfc.fa: %.fa.gz
	#bfc: malloc.c:2839: mremap_chunk: Assertion `((size + offset) & (_rtld_global_ro._dl_pagesize - 1)) == 0' failed.
	#bfc -t$j $< >$@
	bfc $< >$@

# SGA

# Preprocess reads using SGA
%.preprocess.fa: data/%.fa
	sga preprocess -v -o $@ $<

# Index reads using SGA
%.bwt: %.fa
	sga index -v -t$j $<

# Remove subsumed reads using SGA
%.rmdup.fa %.rmdup.bwt: %.fa %.bwt
	sga rmdup -v -t$j $<

# Overlap reads using SGA
sga-m500/%.asqg.gz: %.fa %.bwt
	sga overlap -v -t$j -m500 $<

# Assemble reads using SGA
%.sga-contigs.fa: %.asqg.gz
	sga assemble -v -m500 -o $*.sga $<

# Assemble contigs using Celera Assembler

# Uncompress a FASTQ file
%.fa: %.fa.gz
	gunzip -c $< >$@

# Create a FRG file from a FASTQ file
celera/%.frg: data/%.fa
	mkdir -p $(@D)
	fastqToCA -libraryname nanopore -technology pacbio-corrected -reads $< >$@

# Assemble contigs using Celera Assembler
celera/scerevisiae.scf.fasta: celera/c-elegans_LongRead_500_1499nt.frg celera/c-elegans_LongRead.frg
	runCA -d $(@D) -p scerevisiae $^

# Create a FRG file from a FASTQ file
%.frg: %.fa.gz
	fastqToCA -libraryname nanopore -technology pacbio-corrected -reads $< >$@

# Assemble contigs using Celera Assembler using a spec file
%/scerevisiae.scf.fasta: %.spec data/W303_ONT_Nanocorr_Corrected_reads.frg
	runCA -d $(@D) -p scerevisiae -s $^

# ABySS

# Assemble unitigs using a standard de Bruijn Graph
k$k/%-unitigs.fa:
	mkdir -p k$k
	$(shell which time) -p abyss-pe -C k$k \
		name=scerevisiae \
		k=$k e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		long=nanopore nanopore='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		unitigs

# Assemble scaffolds using a standard de Bruijn Graph
k$k/%-long-scaffs.fa: k$k/%-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k \
		name=scerevisiae \
		k=$k e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		long=nanopore nanopore='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz'

# Assemble using a paired dBG

# Assemble unitigs
k$k-K$K/%-unitigs.fa:
	mkdir -p k$k-K$K
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=scerevisiae \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		long=nanopore nanopore='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		unitigs

# Scaffold with Nanopore data
k$k-K$K/%-long-scaffs.fa: k$k-K$K/%-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=scerevisiae \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		long=nanopore nanopore='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz'

# Create symlinks for the unitigs
k$k-K$K-scaff/scerevisiae-unitigs.fa: k$k-K$K/scerevisiae-unitigs.fa
	mkdir -p $(@D)
	ln -sf ../$(<D)/{scerevisiae-1.fa,scerevisiae-1.dot,scerevisiae-1.path,scerevisiae-2.fa,scerevisiae-2.dot1,scerevisiae-2.dot,scerevisiae-2.path,scerevisiae-3.fa,scerevisiae-3.dot} $(@D)/
	ln -sf scerevisiae-3.fa $(@D)/scerevisiae-unitigs.fa

# Scaffold
k$k-K$K-scaff/%-scaffolds.fa: k$k-K$K-scaff/%-unitigs.fa
	$(shell which time) -p abyss-pe -C $(@D) \
		name=scerevisiae \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		long=nanopore nanopore='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		scaffolds

# Scaffold with Nanopore data
k$k-K$K-scaff/%-long-scaffs.fa: k$k-K$K-scaff/%-scaffolds.fa
	$(shell which time) -p abyss-pe -C $(@D) \
		name=scerevisiae \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz' \
		long=nanopore nanopore='../data/W303_ONT_Nanocorr_Corrected_reads.fa.gz'

# Construct a Bloom filter
scerevisiae.k%.bloom: data/c-elegans_LongRead_500_1499nt.fastq.gz data/c-elegans_LongRead.fastq.gz
	abyss-bloom build -v -k$* -j$j -b500M -l2 $@ $^

# Fill in gaps using ABySS-sealer
k$k-K$K-sealer/%-scaffolds.fa: k$k-K$K-scaff/%-scaffolds.fa \
		%.k25.bloom \
		%.k50.bloom \
		%.k75.bloom \
		%.k100.bloom \
		%.k150.bloom \
		%.k200.bloom \
		%.k250.bloom \
		%.k300.bloom \
		%.k400.bloom \
		%.k500.bloom
	mkdir -p k$k-K$K-sealer
	ln -sf ../k$k-K$K-scaff/$*-unitigs.fa k$k-K$K-sealer/
	abyss-sealer -v -j$j \
		--print-flanks \
		--max-frag=2000 -L512 -k25 -k50 -k75 -k100 -k150 -k200 -k250 -k300 -k400 -k500 \
		-o k$k-K$K-sealer/$* -t k$k-K$K-sealer/$*_trace.tsv -S $< \
		$(addprefix -i , $(wordlist 2, 99, $^)) \
		/dev/null
	ln -sf $*_scaffold.fa k$k-K$K-sealer/$*-scaffolds.fa

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	abyss-fatoagp -f $@ $< >$@.agp

# Assemble the reference genome
$(ref)-k%/scerevisiae-1.fa: $(ref).fa
	mkdir -p $(ref)-k$*
	ABYSS -v -k$* -e0 -t0 -c0 $< -o $@ -s $(ref)-k$*/scerevisiae-bubbles.fa 2>&1 |tee $@.log

# Make assembly symlinks

$(ref)-%/scerevisiae.bwa.fa: $(ref).W303_ONT_Nanocorr_Corrected_reads.sort.bam.%.bed.fa
	mkdir -p $(@D)
	ln -sf ../$< $@

sga-m500/sga-contigs.fa: sga-m500/c-elegans_LongRead_all.bfc.preprocess.rmdup.sga-contigs.fa
	ln -sf $(<F) $@

nanocorr/scerevisiae-polished.fa.gz: data/W303_ONT_Assembly_CA_polished.fa.gz
	mkdir -p $(@D)
	ln -sf ../$< $@

# Calculate assembly statistics
assembly-stats.tsv: \
		$(ref)-k64/scerevisiae-1.fa \
		$(ref)-k128/scerevisiae-1.fa \
		$(ref)-k256/scerevisiae-1.fa \
		$(ref)-k512/scerevisiae-1.fa \
		$(ref)-depth1+/scerevisiae.bwa.fa \
		$(ref)-depth2+/scerevisiae.bwa.fa \
		nanocorr/scerevisiae-polished.fa.gz \
		k512/scerevisiae-unitigs.fa \
		k512-K96/scerevisiae-unitigs.fa \
		k512-K96-scaff/scerevisiae-unitigs.fa \
		k512-K96-scaff/scerevisiae-scaffolds.fa \
		k512/scerevisiae-long-scaffs.fa \
		k512-K96/scerevisiae-long-scaffs.fa \
		k512-K96-scaff/scerevisiae-long-scaffs.fa
	abyss-fac -e12157105 -s1000 $^ >$@

# Analayze the assemblies using QUAST
%.quast/transposed_report.tsv: \
		k512/scerevisiae-long-scaffs.fa \
		k512-K96/scerevisiae-long-scaffs.fa \
		k512-K96-scaff/scerevisiae-long-scaffs.fa
	quast.py -T$j -esL -o $*.quast -R $(ref).fa -G $(ref).gff $^

# Convert TSV to Markdown
%.tsv.md: %.tsv
	abyss-tabtomd $< >$@

# Render HTML from Markdown
%.html: %.md
	pandoc -o $@ $<

# Render HTML from RMarkdown
%.html: %.rmd
	Rscript -e 'rmarkdown::render("$<", "html_document")'

# Dependencies

assembly-stats.html: assembly-stats.tsv

quast-stats.html: scerevisiae.quast/transposed_report.tsv
