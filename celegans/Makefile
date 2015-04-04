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

# Reference genome
ref=WBcel235

# Targets

all: \
	assembly-stats.tsv \
	assembly-stats.html \
	celegans.quast/transposed_report.tsv \
	quast-stats.html

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Download data

# Download the reference genome
$(ref).fa:
	curl ftp://ftp.ensembl.org/pub/release-78/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz |gunzip -c >$@

# Download gene annotations from Ensembl in GTF format
$(ref).gtf:
	curl ftp://ftp.ensembl.org/pub/release-78/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.78.gtf.gz |gunzip -c >$@

# Download gene annotations WormBase in GFF format
$(ref).gff:
	curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.WS240.annotations.gff3.gz |gunzip -c >$@

# Download the Moleculo reads shorter than 1500 bp
data/c-elegans_LongRead_500_1499nt.fastq.gz:
	curl -o $@ https://api.basespace.illumina.com/v1pre3/files/691535736/content

# Download the Moleculo reads at least 1500 bp
data/c-elegans_LongRead.fastq.gz:
	curl -o $@ https://api.basespace.illumina.com/v1pre3/files/691535735/content

# Concatenate the Moleculo data
data/c-elegans_LongRead_all.fastq.gz: data/c-elegans_LongRead_500_1499nt.fastq.gz data/c-elegans_LongRead.fastq.gz
	cat $^ >$@

# ERX046316 2.5 kbp mate-pair
data/ERR068453_%.fastq.bz2:
	curl -o $@ ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/ERA081/ERA081001/ERX046316/$(@F)

# DRX007634 4.5 kbp mate-pair
data/DRR008445_%.fastq.bz2:
	curl -o $@ ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA000/DRA000967/DRX007634/$(@F)

# BWA

# Index the reference
%.fa.bwt: %.fa
	bwa index $<

# Align reads to the reference using BWA MEM
$(ref).%.sam: data/%.fastq.gz $(ref).fa.bwt
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

# Assemble contigs using Celera Assembler

# Create a FRG file from a FASTQ file
celera/%.frg: data/%.fastq.gz
	mkdir -p $(@D)
	fastqToCA -libraryname moleculo -technology moleculo -reads $< >$@

# Assemble contigs using Celera Assembler
celera/celegans.scf.fasta: celera/c-elegans_LongRead_500_1499nt.frg celera/c-elegans_LongRead.frg
	$(shell which time) -p \
	runCA -d $(@D) -p celegans $^

# ABySS

# Merge overlapping reads using abyss-mergepairs
%_merged.fastq: %_1.fq.gz %_2.fq.gz
	abyss-mergepairs -q 15 -v -o $* $^ >$*_merged.tsv

# Assemble unitigs using a standard de Bruijn Graph
k$k/%-unitigs.fa:
	mkdir -p k$k
	$(shell which time) -p abyss-pe -C k$k \
		name=celegans \
		k=$k e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		pe='DRX007634' DRX007634='../data/DRR008445_1.fastq.bz2 ../data/DRR008445_2.fastq.bz2' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz' \
		unitigs

# Assemble contigs using a standard de Bruijn Graph
k$k/%-contigs.fa: k$k/%-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k \
		name=celegans \
		k=$k e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		pe='DRX007634' DRX007634='../data/DRR008445_1.fastq.bz2 ../data/DRR008445_2.fastq.bz2' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz' \
		contigs

# Assemble scaffolds using a standard de Bruijn Graph
k$k/%-scaffolds.fa: k$k/%-contigs.fa
	$(shell which time) -p abyss-pe -C k$k \
		name=celegans \
		k=$k e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		pe='DRX007634' DRX007634='../data/DRR008445_1.fastq.bz2 ../data/DRR008445_2.fastq.bz2' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz' \
		scaffolds

# Assemble scaffolds using a standard de Bruijn Graph
k$k/%-long-scaffs.fa: k$k/%-scaffolds.fa
	$(shell which time) -p abyss-pe -C k$k \
		name=celegans \
		k=$k e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		pe='DRX007634' DRX007634='../data/DRR008445_1.fastq.bz2 ../data/DRR008445_2.fastq.bz2' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz'

# Assemble using a paired dBG

# Assemble unitigs
k$k-K$K/%-unitigs.fa:
	mkdir -p k$k-K$K
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=celegans \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		long=moleculo moleculo='../data/c-elegans_LongRead.fastq.gz' \
		unitigs

# Scaffold with Moleculo data
k$k-K$K/%-long-scaffs.fa: k$k-K$K/%-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=celegans \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz'

# Create symlinks for the unitigs
k$k-K$K-scaff/celegans-unitigs.fa: k$k-K$K/celegans-unitigs.fa
	mkdir -p $(@D)
	ln -sf ../$(<D)/{celegans-1.fa,celegans-1.dot,celegans-1.path,celegans-2.fa,celegans-2.dot1,celegans-2.dot,celegans-2.path,celegans-3.fa,celegans-3.dot} $(@D)/
	ln -sf celegans-3.fa $(@D)/celegans-unitigs.fa

# Scaffold with mate-pair data
k$k-K$K-scaff/%-scaffolds.fa: k$k-K$K-scaff/%-unitigs.fa
	$(shell which time) -p abyss-pe -C $(@D) \
		name=celegans \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		mp='DRX007634' DRX007634='../data/DRR008445_1.fastq.bz2 ../data/DRR008445_2.fastq.bz2' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz' \
		scaffolds

# Scaffold with Moleculo data
k$k-K$K-scaff/%-long-scaffs.fa: k$k-K$K-scaff/%-scaffolds.fa
	$(shell which time) -p abyss-pe -C $(@D) \
		name=celegans \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s n=$n S=$S N=$N j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		mp='DRX007634' DRX007634='../data/DRR008445_1.fastq.bz2 ../data/DRR008445_2.fastq.bz2' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz'

# Construct a Bloom filter
celegans.k%.bloom: data/c-elegans_LongRead_500_1499nt.fastq.gz data/c-elegans_LongRead.fastq.gz
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
		WBcel235-k192/celegans-1.fa \
		WBcel235-k256/celegans-1.fa \
		WBcel235-k512/celegans-1.fa \
		bwa-depth1+/WBcel235.bwa.fa \
		bwa-depth2+/WBcel235.bwa.fa \
		kraghavan-k64/celegans-unitigs.fa \
		kraghavan-k64/celegans-contigs.fa \
		kraghavan-k64/celegans-scaffolds.fa \
		k192/celegans-unitigs.fa \
		k192/celegans-contigs.fa \
		k192/celegans-scaffolds.fa \
		k192/celegans-long-scaffs.fa \
		k192/celegans-scaftigs.fa \
		k512/celegans-unitigs.fa \
		k512/celegans-contigs.fa \
		k512/celegans-scaffolds.fa \
		k512/celegans-long-scaffs.fa \
		k512/celegans-scaftigs.fa \
		k512-K96/celegans-unitigs.fa \
		k512-K96/celegans-contigs.fa \
		k512-K96/celegans-scaffolds.fa \
		k512-K96/celegans-long-scaffs.fa \
		k512-K96-scaff/celegans-unitigs.fa \
		k512-K96-scaff/celegans-scaffolds.fa \
		k512-K96-scaff/celegans-long-scaffs.fa \
		k512-K96-sealer/celegans-unitigs.fa \
		k512-K96-sealer/celegans-scaffolds.fa \
		k512-K96-sealer/celegans-scaftigs.fa
	abyss-fac -e100286401 -s1000 $^ >$@

# Analayze the assemblies using QUAST
%.quast/transposed_report.tsv: \
		kraghavan-k64/celegans-scaffolds.fa \
		k192/celegans-long-scaffs.fa \
		k512/celegans-long-scaffs.fa \
		k512-K96/celegans-long-scaffs.fa \
		k512-K96-scaff/celegans-long-scaffs.fa \
		k512-K96-sealer/celegans-scaffolds.fa
	quast.py -T$j -esL -o $*.quast -R $(ref).fa -G $(ref).gff $^

# Convert TSV to Markdown
%.tsv.md: %.tsv
	abyss-tabtomd $< >$@

# Render HTML and Markdown from RMarkdown
%.html %.md: %.rmd
	Rscript -e 'rmarkdown::render("$<", "html_document")'

# Render HTML from Markdown
%.html: %.md

# Change the resolution of the figure to 300 DPI
%.png: %_files/figure-html/plot-NG50-1.png
	convert $< -units PixelsPerInch -density 300 $@
	pandoc -o $@ $<

# Dependencies

assembly-stats.html assembly-stats.md: assembly-stats.tsv

quast-stats.html quast-stats.md: scerevisiae.quast/transposed_report.tsv
