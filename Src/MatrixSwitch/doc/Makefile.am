#
# Makefile for the documentation of MatrixSwitch
#

# Markdown sources
msw_mddocs = \
  configuring-MatrixSwitch.md \
  contributing-to-MatrixSwitch.md \
  hacking-the-build-system.md

# HTML documents
msw_htmldocs = \
  configuring-MatrixSwitch.html \
  contributing-to-MatrixSwitch.html \
  hacking-the-build-system.html

                    ########################################

# User-defined rules
MARKDOWN = $(top_srcdir)/config/scripts/markdown
MDFLAGS  = 

SUFFIXES = .md .html
.md.html:
	sed -e "s/@TITLE@/$$(head -n 1 $<)/" \
	  $(srcdir)/data/msw_header.html >$@
	$(MARKDOWN) $(MDFLAGS) $< >>$@
	cat $(srcdir)/data/msw_footer.html >>$@

# Targets
htmldocdir = $(docdir)/html
htmldoc_DATA = $(msw_htmldocs)

# Inform Automake
CLEANFILES = $(msw_htmldocs)
EXTRA_DIST = data/msw_header.html data/msw_footer.html $(msw_mddocs)
