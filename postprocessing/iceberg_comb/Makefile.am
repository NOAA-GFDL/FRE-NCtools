bin_SCRIPTS = iceberg_comb.sh

# Documentation settings
A2X = a2x
ASCIIDOC_CFG = $(srcdir)/asciidoc.conf

if WITH_ASCIIDOC
man_MANS = iceberg_comb.sh.1
else
man_MANS =
endif

.txt.1:
	$(A2X) -v -L\
               -D $(builddir) \
               -d manpage \
               -f manpage \
               -afre_version=$(PACKAGE_VERSION) \
               --asciidoc-opts="-f $(ASCIIDOC_CFG)" \
               $<

clean-local:
	$(RM) *.1
