EPStoPDF = epstopdf

ALL: paper.ps paper.pdf 

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

paper.dvi: paper.tex references.bib $(eps_source)
	latex paper.tex < /dev/null
	bibtex paper < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null

paper.pdf: paper.tex references.bib $(pdf_source)
	pdflatex paper.tex < /dev/null
	bibtex paper < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null

pdf:	paper.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) paper.dvi 
	$(RM) paper.blg paper.log
	$(RM) paper.aux paper.ps paper.bbl
	$(RM) *~

.PHONY: clean
