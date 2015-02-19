# Number of threads
j=12

# ABySS parameters
k=512
K=96
l=40
s=500
S=500
e=2
E=0
c=1.01

# Reference genome
ref=WBcel235

# Targets

all: \
	WBcel235.c-elegans_LongRead_all.sort.bam.bai \
	k$k-K$K/celegans-unitigs.fa \
	k$k-K$K/celegans-long-scaffs.fa \
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

# Scaffold with Moleculo data
k$k-K$K/%-long-scaffs.fa: k$k-K$K/%-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=celegans \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s S=$S j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz'

# Create symlinks for the unitigs
k$k-K$K-scaff/celegans-unitigs.fa: k$k-K$K/celegans-unitigs.fa
	mkdir -p $(@D)
	ln -sf ../$(<D)/{celegans-1.fa,celegans-1.dot,celegans-1.path,celegans-2.fa,celegans-2.dot1,celegans-2.dot,celegans-2.path,celegans-3.fa,celegans-3.dot} $(@D)/
	ln -sf celegans-3.fa $(@D)/celegans-unitigs.fa

# Scaffold with mate-pair and Moleculo data
k$k-K$K-scaff/%-long-scaff.fa: k$k-K$K-scaff/%-unitigs.fa
	$(shell which time) -p abyss-pe -C $(@D) \
		name=celegans \
		k=$k K=$K e=$e E=$E c=$c l=$l s=$s S=$S j=$j v=-v \
		se='../data/c-elegans_LongRead_500_1499nt.fastq.gz ../data/c-elegans_LongRead.fastq.gz' \
		mp='ERX046316 DRX007634' \
		ERX046316='../data/ERR068453_1.fastq.bz2 ../data/ERR068453_2.fastq.bz2' \
		DRX007634='../data/DRR008445_1.fastq.bz2 ../data/DRR008445_2.fastq.bz2' \
		long=moleculo moleculo='../data/c-elegans_LongRead_all.fastq.gz'

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
		k512-K96-scaff/celegans-unitigs.fa \
		kraghavan-k64/celegans-contigs.fa \
		kraghavan-k64/celegans-scaffolds.fa \
		k512-K96-scaff/celegans-scaffolds.fa \
		k512/celegans-long-scaffs.fa \
		k512-K96/celegans-long-scaffs.fa \
		k512-K96-scaff/celegans-long-scaffs.fa
	abyss-fac -e100286401 -s1000 $^ >$@

# Analayze the assemblies using QUAST
%.quast/transposed_report.tsv: \
		kraghavan-k64/celegans-scaffolds.fa \
		k512/celegans-long-scaffs.fa \
		k512-K96/celegans-long-scaffs.fa \
		k512-K96-scaff/celegans-long-scaffs.fa
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

quast-stats.html: celegans.quast/transposed_report.tsv
