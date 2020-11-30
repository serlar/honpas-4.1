#
# Makefile for hybrid XC diagrams
#

# Diagram utilities
GV_DOT = dot
GV_DOTFLAGS =

# Diagram building rules
.SUFFIXES:
.SUFFIXES: .pdf .svg .dot

.dot.pdf:
	$(GV_DOT) $(GV_DOTFLAGS) -T pdf -o $@ $<

.dot.svg:
	$(GV_DOT) $(GV_DOTFLAGS) -T svg -o $@ $<

# Diagrams to build
diagrams = \
  siesta-honpas-files-32p5-jul17-nov17.pdf \
  siesta-honpas-files-32p5-jul17-nov17.svg \
  siesta-honpas-files-41b3-jul17-nov17.pdf \
  siesta-honpas-files-41b3-jul17-nov17.svg

# Diagram targets
all all_targets: $(diagrams)

clean distclean:
	rm -f $(diagrams)

