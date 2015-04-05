all: assembly-stats.html

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:

%.html: %.rmd
	Rscript -e 'rmarkdown::render("$<", "html_document")'
