# Makefile for creating the Man page

VERSION = "Bronx"
DATE := $(shell date +"%Y %B %d")

A2X = a2x
ASCIIDOC_CFG = asciidoc.conf

TXTS = iceberg_comb.txt
TXT_MAN1 = $(patsubst %.txt,%.1,$(TXTS))

all: txtman

txtman: $(TXT_MAN1)

clean:
	$(RM) *.[1-9]
	$(RM) *.html
	$(RM) *.tmp
	$(RM) *.tex
	$(RM) *.pdf
	$(RM) *.aux
	$(RM) *.idx
	$(RM) *.log
	$(RM) *.toc
	$(RM) *.xml

%.1: %.txt $(ASCIIDOC_CFG)
	$(A2X) -v \
               -d manpage \
               -f manpage \
               -afre_version=$(VERSION) \
               --asciidoc-opts="-f $(ASCIIDOC_CFG)" \
               $<
