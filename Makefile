# ABySS parameters
k=364
K=32

# Reference genome
ref=NC_000913

# Targets

all: \
	ecoli_merged.fastq \
	k$k/ecoli-scaffolds.fa \
	k$k/ecoli-scaftigs.fa \
	k$k-K$K/ecoli-unitigs.fa \
	k$k-K$K/pe600-3.dist \
	k$k-K$K-scaff/ecoli-scaffolds.fa \
	k$k-K$K-sealer/ecoli-scaffolds.fa \
	k$k-K$K-sealer/ecoli-scaftigs.fa \
	ecoli-assembly-stats.tsv \
	ecoli-assembly-stats.html

other: \
	k$k-K$K/ecoli-contigs.fa \
	k$k-K$K/ecoli-scaffolds.fa \

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Merge overlapping reads using abyss-mergepairs
%_merged.fastq: %_1.fq.gz %_2.fq.gz
	abyss-mergepairs -q 15 -v -o $* $^ >$*_merged.tsv

# Assemble scaffolds using a standard de Bruijn Graph
k$k/%-scaffolds.fa: %_merged.fastq
	mkdir -p k$k
	$(shell which time) -p abyss-pe -C k$k \
		name=ecoli \
		k=$k l=40 s=1000 v=-v \
		se='../ecoli_merged.fastq ../ecoli_reads_1.fastq ../ecoli_reads_2.fastq' \
		pe='pe600' \
		pe600='../ecoli_1.fq.gz ../ecoli_2.fq.gz'

# Assemble using a paired dBG

# Assemble unitigs
k$k-K$K/%-unitigs.fa: %_merged.fastq
	mkdir -p k$k-K$K
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=ecoli \
		k=$k K=$K l=40 s=1000 v=-v \
		se='../ecoli_merged.fastq ../ecoli_reads_1.fastq ../ecoli_reads_2.fastq' \
		pe='pe600' \
		pe600='../ecoli_1.fq.gz ../ecoli_2.fq.gz' \
		unitigs

# Estimate distances
k$k-K$K/pe600-3.dist: k$k-K$K/ecoli-unitigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=ecoli \
		k=$k K=$K l=40 s=1000 v=-v \
		se='../ecoli_merged.fastq ../ecoli_reads_1.fastq ../ecoli_reads_2.fastq' \
		pe='pe600' \
		pe600='../ecoli_1.fq.gz ../ecoli_2.fq.gz' \
		pe-sam pe600-3.dist

# Assemble contigs
k$k-K$K/ecoli-contigs.fa: k$k-K$K/pe600-3.dist
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=ecoli \
		k=$k K=$K l=40 s=1000 v=-v \
		se='../ecoli_merged.fastq ../ecoli_reads_1.fastq ../ecoli_reads_2.fastq' \
		pe='pe600' \
		pe600='../ecoli_1.fq.gz ../ecoli_2.fq.gz' \
		contigs

# Assemble scaffolds
k$k-K$K/ecoli-scaffolds.fa: k$k-K$K/ecoli-contigs.fa
	$(shell which time) -p abyss-pe -C k$k-K$K \
		name=ecoli \
		k=$k K=$K l=40 s=1000 v=-v \
		se='../ecoli_merged.fastq ../ecoli_reads_1.fastq ../ecoli_reads_2.fastq' \
		pe='pe600' mp='pe600' \
		pe600='../ecoli_1.fq.gz ../ecoli_2.fq.gz' \
		scaffolds scaftigs stats

# Scaffold without first contigging
k$k-K$K-scaff/ecoli-scaffolds.fa: k$k-K$K/pe600-3.dist
	mkdir -p k$k-K$K-scaff
	ln -sf ../k$k-K$K/ecoli-3.fa k$k-K$K-scaff/
	ln -sf ecoli-3.fa k$k-K$K-scaff/ecoli-unitigs.fa
	ln -sf ecoli-3.fa k$k-K$K-scaff/ecoli-6.fa
	ln -sf ../k$k-K$K/ecoli-3.dot k$k-K$K-scaff/ecoli-6.dot
	ln -sf ../k$k-K$K/pe600-3.dist k$k-K$K-scaff/pe600-6.dist
	abyss-todot -e k$k-K$K-scaff/ecoli-6.fa k$k-K$K-scaff/pe600-6.dist >k$k-K$K-scaff/pe600-6.dist.dot
	$(shell which time) -p abyss-pe -C k$k-K$K-scaff \
		name=ecoli \
		k=$k K=$K l=40 s=1000 v=-v \
		se='../ecoli_merged.fastq ../ecoli_reads_1.fastq ../ecoli_reads_2.fastq' \
		pe='pe600' \
		pe600='../ecoli_1.fq.gz ../ecoli_2.fq.gz' \
		-o ecoli-6.fa -o ecoli-6.dot -o pe600-6.dist.dot \
		scaffolds scaftigs

# Scaffold using abyss-scaffold
foo:
	abyss-scaffold -v   -k364 -s1000 -n10 -g ecoli-6.path.dot  ecoli-6.dot pe600-6.dist.dot >ecoli-6.path
	PathConsensus -v --dot -k364  -p0.9 -s ecoli-7.fa -g ecoli-7.dot -o ecoli-7.path ecoli-6.fa ecoli-6.dot ecoli-6.path
	cat ecoli-6.fa ecoli-7.fa \
		|MergeContigs -v  -k364 -o ecoli-8.fa - ecoli-7.dot ecoli-7.path
	ln -sf ecoli-8.fa ecoli-scaffolds.fa

# Fill in gaps using ABySS-sealer
k$k-K$K-sealer/%-scaffolds.fa: k$k-K$K-scaff/%-scaffolds.fa
	mkdir -p k$k-K$K-sealer
	ln -sf ../k$k-K$K-scaff/$*-unitigs.fa k$k-K$K-sealer/
	abyss-sealer -v -j12 \
		--print-flanks \
		-L364 -k50 -k100 -k150 -k200 -k250 -k300 -k350 -k364 \
		-o k$k-K$K-sealer/$* -S $< \
		$*_merged.fastq $*_reads_1.fastq $*_reads_2.fastq
	ln -s $*_scaffold.fa k$k-K$K-sealer/$*-scaffolds.fa

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	abyss-fatoagp -f $@ $< >$@.agp

# Assemble the reference genome
$(ref)-k%/ecoli-1.fa: $(ref).fa
	mkdir -p $(ref)-k$*
	ABYSS -v -k$* -e0 -t0 -c0 $< -o $@ -s $(ref)-k$*/ecoli-bubbles.fa 2>&1 |tee $@.log

# Calculate assembly statistics
%-assembly-stats.tsv: k$k-K$K-sealer/%-scaffolds.fa k$k-K$K-sealer/%-scaftigs.fa
	abyss-fac \
		$(ref)-k64/ecoli-1.fa \
		$(ref)-k364/ecoli-1.fa \
		k$k/ecoli-unitigs.fa \
		k$k/ecoli-contigs.fa \
		k$k/ecoli-scaffolds.fa \
		k$k/ecoli-scaftigs.fa \
		k416/ecoli-unitigs.fa \
		k416/ecoli-contigs.fa \
		k416/ecoli-scaffolds.fa \
		k416/ecoli-scaftigs.fa \
		k$k-K$K/ecoli-unitigs.fa \
		k$k-K$K/ecoli-contigs.fa \
		k$k-K$K/ecoli-scaffolds.fa \
		k$k-K$K/ecoli-scaftigs.fa \
		k$k-K$K-sealer/ecoli-unitigs.fa \
		$^ \
		>$@

# Convert TSV to Markdown
%.tsv.md: %.tsv
	abyss-tabtomd $< >$@

# Render HTML from Markdown
%.html: %.md
	pandoc -o $@ $<

# Render HTML from RMarkdown
%.html: %.rmd %.tsv
	Rscript -e 'rmarkdown::render("$<", "html_document")'
