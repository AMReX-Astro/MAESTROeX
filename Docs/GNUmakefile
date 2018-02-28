EPStoPDF = epstopdf

ALL: MAESTROeXUsersGuide.pdf 

eps_source = $(wildcard *.eps) \
             $(wildcard lo_density/*.eps) \
             $(wildcard flowchart/*.eps) \
             $(wildcard architecture/*.eps) \
             $(wildcard unit_tests/*.eps) \
             $(wildcard initial_models/*.eps) \
             $(wildcard rotation/*.eps) \
             $(wildcard radial_base/*.eps) \
             $(wildcard spherical_basestate/*.eps) \
             $(wildcard visualization/*.eps) \
             $(wildcard planar_invsq_base/*.eps) \
             $(wildcard pert_notes/*.eps) \
             $(wildcard SDC/*.eps) \
             $(wildcard getting_started/*.eps) \
             $(wildcard mg/*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = maestro_doc.bib \
                  introduction/introduction.tex \
                  preface/preface.tex \
                  flowchart/flowchart.tex \
                  getting_started/getting_started.tex \
                  unit_tests/unit_tests.tex \
                  architecture/architecture.tex \
                  runtime_parameters/runtime_parameters.tex \
                  volume_discrepancy/volume_discrepancy.tex \
                  eos_notes/eos_notes.tex \
                  networks/networks.tex \
                  eta_notes/eta.tex \
                  pert_notes/pert.tex \
                  Godunov_notes/Godunov.tex \
                  initial_models/initial_models.tex \
                  lo_density/lo_density.tex \
                  rotation/rotation.tex \
                  radial_base/radial_base.tex \
                  spherical_basestate/basestate.tex \
                  FAQ/faq.tex \
                  analysis/analysis.tex \
                  thermo_notes/thermo_notes.tex \
                  visualization/visualization.tex \
                  SDC/sdc.tex \
                  makefiles/makefiles.tex \
                  managing_jobs/managingjobs.tex \
                  mg/mg.tex \
                  MESA/instructions.tex


MAESTROeXUsersGuide.pdf: MAESTROeXUsersGuide.tex maestrosymbols.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex MAESTROeXUsersGuide.tex < /dev/null
	bibtex MAESTROeXUsersGuide
	pdflatex MAESTROeXUsersGuide.tex < /dev/null
	makeindex MAESTROeXUsersGuide < /dev/null
	pdflatex MAESTROeXUsersGuide.tex < /dev/null
	pdflatex MAESTROeXUsersGuide.tex < /dev/null

pdf:	MAESTROeXUsersGuide.pdf 


%.pdf: %.eps
	$(EPStoPDF) $<


clean:
	$(RM) $(pdf_source) MAESTROeXUsersGuide.dvi 
	$(RM) MAESTROeXUsersGuide.blg MAESTROeXUsersGuide.log 
	$(RM) MAESTROeXUsersGuide.toc MAESTROeXUsersGuide.lof
	$(RM) MAESTROeXUsersGuide.idx  MAESTROeXUsersGuide.ind  
	$(RM) MAESTROeXUsersGuide.ilg  MAESTROeXUsersGuide.lot
	$(RM) MAESTROeXUsersGuide.aux MAESTROeXUsersGuide.ps MAESTROeXUsersGuide.bbl
	$(RM) *~

.PHONY: clean
