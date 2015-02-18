# Number of threads
j=12

# ABySS parameters
k=512
K=96
l=40
s=500
S=500
e=2
E=1
c=2

# Reference genome
ref=WBcel235

# Targets

all: \
	WBcel235.c-elegans_LongRead_all.sort.bam.bai \
	k$k-K$K/celegans-unitigs.fa \
	k$k-K$K/celegans-long-scaffs.fa \
	assembly-stats.tsv \
	assembly-stats.html

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Download data

data/ERR068453_%.fastq.bz2:
	curl -o $@ ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/ERA081/ERA081001/ERX046316/$(@F)

data/DRR008445_%.fastq.bz2:
	curl -o $@ ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA000/DRA000967/DRX007634/$(@F)

# BWA

# Index the reference
%.fa.bwt: %.fa
	bwa index $<

# Align reads to the reference using BWA MEM
$(ref).%.sam: %.sam $(ref).fa.bwt
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
%.bfc.fastq: %.fastq.gz
	#bfc: malloc.c:2839: mremap_chunk: Assertion `((size + offset) & (_rtld_global_ro._dl_pagesize - 1)) == 0' failed.
	#bfc -t$j $< >$@
	bfc $< >$@

# SGA

# Preprocess reads using SGA
%.preprocess.fastq: data/%.fastq
	sga preprocess -v -o $@ $<

# Index reads using SGA
%.bwt: %.fastq
	sga index -v -t$j $<

# Remove subsumed reads using SGA
%.rmdup.fa %.rmdup.bwt: %.fastq %.bwt
	sga rmdup -v -t$j $<

# Overlap reads using SGA
sga-m500/%.asqg.gz: %.fa %.bwt
	sga overlap -v -t$j -m500 $<

# Assemble reads using SGA
%.sga-contigs.fa: %.asqg.gz
	sga assemble -v -m500 -o $*.sga $<

# Merge overlapping reads using abyss-mergepairs
%_merged.fastq: %_1.fq.gz %_2.fq.gz
	abyss-mergepairs -q 15 -v -o $* $^ >$*_merged.tsv

# Assemble unitigs using a standard de Bruijn Graph
k$k/%-unitigs.fa:
	mkdir -p k$k
	$(shell which time) -p abyss-pe -C k$k \
		name=celegans \
		k=$k e=$e E=$E c=$c l=$l s=$s S=$S j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz' \
		unitigs

# Assemble scaffolds using a standard de Bruijn Graph
k$k/%-long-scaffs.fa: k$k/%-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k \
		name=celegans \
		k=$k e=$e E=$E c=$c l=$l s=$s S=$S j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz'

# Assemble using a paired dBG

# Assemble unitigs
k$k-K$K/%-unitigs.fa:
	mkdir -p k$k-K$K
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=celegans \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s S=$S j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		long=moleculo moleculo='../data/c-elegans_LongRead.fastq.gz' \
		unitigs

# Assemble scaffolds
k$k-K$K/%-long-scaffs.fa: k$k-K$K/%-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=celegans \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s S=$S j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz'

# Estimate distances
k$k-K$K/moleculo-3.dist: k$k-K$K/celegans-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=celegans \
		k=$k K=$K l=$l s=$s v=-v \
		se='../celegans_merged.fastq ../celegans_reads_1.fastq ../celegans_reads_2.fastq' \
		long='moleculo' \
		moleculo='../celegans_1.fq.gz ../celegans_2.fq.gz' \
		pe-sam moleculo-3.dist

# Assemble contigs
k$k-K$K/celegans-contigs.fa: k$k-K$K/moleculo-3.dist
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=celegans \
		k=$k K=$K l=$l s=$s v=-v \
		se='../celegans_merged.fastq ../celegans_reads_1.fastq ../celegans_reads_2.fastq' \
		long='moleculo' \
		moleculo='../celegans_1.fq.gz ../celegans_2.fq.gz' \
		contigs

# Assemble scaffolds
k$k-K$K/celegans-scaffolds.fa: k$k-K$K/celegans-contigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=celegans \
		k=$k K=$K l=$l s=$s v=-v \
		se='../celegans_merged.fastq ../celegans_reads_1.fastq ../celegans_reads_2.fastq' \
		long='moleculo' \
		moleculo='../celegans_1.fq.gz ../celegans_2.fq.gz' \
		scaffolds scaftigs stats

# Estimate distances between contigs
k$k-K$K-scaff/celegans-4.dot: k$k-K$K/celegans-unitigs.fa
	mkdir -p k$k-K$K-scaff
	ln -sf ../k$k-K$K/celegans-3.fa ../k$k-K$K/celegans-3.dot k$k-K$K-scaff/
	ln -sf celegans-3.fa k$k-K$K-scaff/celegans-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K-scaff \
		name=celegans \
		k=$k K=$K l=$l s=$s v=-v \
		OVERLAP_OPTIONS=--no-scaffold \
		se='../celegans_merged.fastq ../celegans_reads_1.fastq ../celegans_reads_2.fastq' \
		long='moleculo' \
		moleculo='../celegans_1.fq.gz ../celegans_2.fq.gz' \
		-o celegans-1.fa -o celegans-1.dot -o celegans-1.path -o celegans-2.dot1 -o celegans-2.dot -o celegans-3.fa -o celegans-3.dot \
		pe-sam celegans-4.dot

# Scaffold without first contigging
k$k-K$K-scaff/celegans-scaffolds.fa: k$k-K$K-scaff/celegans-4.dot
	ln -sf celegans-3.fa k$k-K$K-scaff/celegans-6.fa
	ln -sf celegans-6.fa k$k-K$K-scaff/celegans-contigs.fa
	ln -sf celegans-4.dot k$k-K$K-scaff/celegans-6.dot
	abyss-todot -e k$k-K$K-scaff/celegans-3.fa k$k-K$K-scaff/moleculo-3.dist >k$k-K$K-scaff/moleculo-3.dist.dot
	ln -sf moleculo-3.dist.dot k$k-K$K-scaff/moleculo-8.dist.dot
	$(shell which time) -p abyss-pe -C k$k-K$K-scaff \
		name=celegans \
		k=$k K=$K l=$l s=$s v=-v \
		se='../celegans_merged.fastq ../celegans_reads_1.fastq ../celegans_reads_2.fastq' \
		long='moleculo' \
		moleculo='../celegans_1.fq.gz ../celegans_2.fq.gz' \
		-o celegans-6.fa -o celegans-6.dot -o moleculo-8.dist.dot \
		scaffolds scaftigs

# Construct a Bloom filter
celegans.k%.bloom: celegans_merged.fastq celegans_reads_1.fastq celegans_reads_2.fastq
	abyss-bloom build -v -k$* -j$j -b500M -l2 $@ $^

# Fill in gaps using ABySS-sealer
k$k-K$K-sealer/%-scaffolds.fa: k$k-K$K-scaff/%-scaffolds.fa \
		%.k50.bloom \
		%.k75.bloom \
		%.k100.bloom \
		%.k150.bloom \
		%.k250.bloom
	mkdir -p k$k-K$K-sealer
	ln -sf ../k$k-K$K-scaff/$*-unitigs.fa k$k-K$K-sealer/
	abyss-sealer -v -j$j \
		--print-flanks \
		--max-frag=2000 -L512 -k50 -k75 -k100 -k150 -k250 \
		-o k$k-K$K-sealer/$* -t k$k-K$K-sealer/$*_trace.tsv -S $< \
		$(addprefix -i , $(wordlist 2, 99, $^)) \
		$*_merged.fastq $*_reads_1.fastq $*_reads_2.fastq
	ln -sf $*_scaffold.fa k$k-K$K-sealer/$*-scaffolds.fa

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	abyss-fatoagp -f $@ $< >$@.agp

# Download the reference genome
$(ref).fa:
	curl ftp://ftp.ensembl.org/pub/release-78/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz |gunzip -c >$@

# Assemble the reference genome
$(ref)-k%/celegans-1.fa: $(ref).fa
	mkdir -p $(ref)-k$*
	ABYSS -v -k$* -e0 -t0 -c0 $< -o $@ -s $(ref)-k$*/celegans-bubbles.fa 2>&1 |tee $@.log

# Make assembly symlinks

bwa-%/WBcel235.bwa.fa: WBcel235.c-elegans_LongRead_all.sort.bam.%.bed.fa
	mkdir -p $(@D)
	ln -sf ../$< $@

sga-m500/sga-contigs.fa: sga-m500/c-elegans_LongRead_all.bfc.preprocess.rmdup.sga-contigs.fa
	ln -sf $(<F) $@

# Calculate assembly statistics
assembly-stats.tsv: \
		WBcel235-k64/celegans-1.fa \
		WBcel235-k128/celegans-1.fa \
		WBcel235-k256/celegans-1.fa \
		WBcel235-k512/celegans-1.fa \
		bwa-depth1+/WBcel235.bwa.fa \
		bwa-depth2+/WBcel235.bwa.fa \
		kraghavan-k64/celegans-unitigs.fa \
		k512/celegans-unitigs.fa \
		k512-K96/celegans-unitigs.fa \
		sga-m500/sga-contigs.fa \
		kraghavan-k64/celegans-contigs.fa \
		kraghavan-k64/celegans-scaffolds.fa \
		k512/celegans-long-scaffs.fa \
		k512-K96/celegans-long-scaffs.fa
	abyss-fac -e100286401 -s1000 $^ >$@

# Analayze the assemblies using QUAST
%.quast/report.tsv: \
		$(ref)-k64/%-1.fa \
		$(ref)-k364/%-1.fa \
		k$k/%-scaffolds.fa \
		k416/%-scaffolds.fa \
		k$k-K$K/%-scaffolds.fa \
		k$k-K$K-scaff/celegans-scaffolds.fa \
		k$k-K$K-sealer/%-scaffolds.fa
	quast.py -fsL -o $*.quast -R $(ref).fa -G $(ref).gff $^

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
