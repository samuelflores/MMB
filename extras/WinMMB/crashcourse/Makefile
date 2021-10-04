MAINFILE = crashcourse

OBJECTS = $(shell echo *.tex)

OBJECTS_TEST = $(addsuffix .t, $(basename $(OBJECTS)))
TESTS = $(addprefix make/, $(OBJECTS_TEST))
TEMP2 := $(shell mkdir make 2>/dev/null)

.PHONY: all
all: $(MAINFILE).pdf

$(MAINFILE).pdf: $(TESTS) $(EXTRA_FILES)
	xelatex $(MAINFILE)
	xelatex $(MAINFILE)

make/%.t: %.tex
	touch $@

.PHONY: clean
clean:
	-rm -f *.aux
	-rm -f *.log
	-rm -f *.toc
	-rm -f *.bbl
	-rm -f *.blg
	-rm -f *.out
	-rm -f make/bib

.PHONY: cleanall
cleanall: clean
	-rm -f *.pdf
	-rm -rf ./make

