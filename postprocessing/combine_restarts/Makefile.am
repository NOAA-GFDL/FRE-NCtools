bin_SCRIPTS = combine_restarts

# Documentation settings
#man_MANS = combine_restarts.1

A2X = a2x
ASCIIDOC_CFG = $(srcdir)/asciidoc.conf

docs-local: combine_restarts.1

install-docs-local: docs-local
	$(INSTALL) -m 755 -d $(DESTDIR)$(mandir)/1
	$(INSTALL_DATA) $(builddir)/combine_restarts.1 $(mandir)/1

combine_restarts.1: $(srcdir)/combine_restarts.txt $(ASCIIDOC_CFG)
	$(A2X) -v \
               -D $(builddir) \
               -d manpage \
               -f manpage \
               -afre_version=$(PACKAGE_VERSION) \
               --asciidoc-opts="-f $(ASCIIDOC_CFG)" \
               $(srcdir)/combine_restarts.txt

clean-local:
	$(RM) *.1
