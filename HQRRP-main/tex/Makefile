#!/usr/bin/env make

LC   = lualatex
BC   = bibtex
LOPT = 
BOPT = 

texname  = PEARC21-CPY21

.PHONY : all bib pdf clean help

all: 
	$(LC) $(LOPT) $(texname)
	$(BC) $(BOPT) $(texname)
	$(LC) $(LOPT) $(texname)
	$(LC) $(LOPT) $(texname)

pdf: 
	$(LC) $(LOPT) $(texname)

bib:
	$(BC) $(BOPT) $(tex1name)

clean:
	rm -f $(texname).{aux,log,bbl,blg,pdf}

help:
	@echo "COMMAND: DESCRIPTION"
	@echo "    all: default -- generates pdf and bib with appr. num. of compiles"
	@echo "    pdf:         -- run latex compiler once on tex"
	@echo "    bib:         -- run bib compiler once"
	@echo "  clean:         -- remove all latex compile related nontex files"
