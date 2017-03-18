# Makefile for creating the Man page

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin

VPATH = $(SRCDIR):.

VERSION = "Bronx"
DATE := $(shell date +"%Y %B %d")

A2X = a2x
ASCIIDOC_CFG = $(SRCDIR)/asciidoc.conf

TARGETS_SHARE = iceberg_comb.sh

TXTS = iceberg_comb.txt
TXT_MAN1 = $(patsubst %.txt,%.1,$(TXTS))

all: 

install: all $(TARGET_SHARE)
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_SHARE)
	install -m 755 -d $(PREFIX_SHARE_EXEC)
	install -m 755 $(TARGETS_SHARE:%=$(SRCDIR)/%) $(PREFIX_SHARE_EXEC)

docs: txtman

install-docs: docs
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_DOCS)
	install -m 755 -d $(PREFIX_DOCS)/man
	install -m 755 -d $(PREFIX_DOCS)/man/1
	install -m 644 $(TXT_MAN1) $(PREFIX_DOCS)/man/1

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
