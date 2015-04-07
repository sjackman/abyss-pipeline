all: assembly-stats.html poster.html

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:

%.html: %.rmd
	Rscript -e 'rmarkdown::render("$<", "html_document")'
