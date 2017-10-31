# Contributing to FRE NCtools

To add a new tool to the FRE NCtools suite of tools, create a
directory in either the [postprocessing](postprocessing) or
[tools](tools) directoy.  In the new directory, include the source
code or scripts for the tool.  Also include a `Makefile` that can be
used to build the package.  Note, the
[fre-nctools-make-package](fre-nctools-make-package) script uses a
separate build directory.  The `Makefile` included must use the
`VPATH` macro to allow make to find the source files.  The `Makefile`
must include the `install`, `docs`, and `install-docs` targets, and
should use the `PREFIX` and `SITE` macros (which are passed in from
the master `Makefile`) to install the compiled binaries into the
`$(PREFIX)/$(SITE)/bin` directory, and the non-compiled scripts into
the `$(PREFIX)/share/bin` directory.  The following `Makefile`
template can be used as a starting point:

```make
# Default compilers
CC       := icc
FC       := ifort

# Bring in setting for the build
# Compilers may be overridden
include ../../build.mk

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin

VPATH = $(SRCDIR):.

# List of targets to build:
TARGETS := <myTargets>
DOC_TARGETS := <docTargets>

# Default build target
all: $(TARGETS)

install: $(TARGETS)
        install -m 755 -d $(PREFIX)
        install -m 755 -d $(PREFIX_EXEC)
        install -m 755 $(TARGETS) $(PREFIX_EXEC)

docs: $(DOC_TARGETS)

install-docs: docs
        install -m 755 -d $(PREFIX)
        install -m 755 -d $(PREFIX_DOCS)
        install -m 755 -d $(PREFIX_DOCS)/man
        install -m 755 -d $(PREFIX_DOCS)/man/1
        install -m 644 $(DOC_TARGETS) $(PREFIX_DOCS)/man/1
```

The master [Makefile](Makefile) must also be updated to include four
new targets: `<tool>:`, `<tool>-docs:`, `<tool>-install:` and
`<tool>-install-docs:`.  The
[fre-nctools-make-package](fre-nctools-make-package) script must also
be updated to include the `<tool>` name in either the `toolsSRC` or
`postpSRC` list.

# Submit a Pull Request / Merge Request

Once the new tool is added to the FRE NCtools repository, submit a
Pull Request (on github) or a Merge Request (if on the GFDL Gitlab
repository).  A member of the Modeling Systems group will review the
request, and will indicate if the package is accepted or if any
changes are needed before being accepted.

# Documentation

Many of the tools in FRE NCtools are still missing documentation.  We
are accepting markdown formatted documentation for any of the FRE
NCtools.  To contribute documentation, please add the markdown version
of the documentation to the tool's directory.  If needed add to the
`Makefile` the build instructions, and submit a Pull/Merge request.
