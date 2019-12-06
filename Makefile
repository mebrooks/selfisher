R=R
# -> you can do    R=R-devel  make ....

PACKAGE=selfisher
# get VERSION from selfisher/DESCRIPTION  
## ("::" = expand only  once, but doesn't work in make <= 3.81)
VERSION := $(shell sed -n '/^Version: /s///p' selfisher/DESCRIPTION)

TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

CPP_SRC := $(PACKAGE)/src/*.cpp

all:
	make enum-update
	make doc-update
	make build-package
	make install
	make pdf

enum-update:: $(PACKAGE)/R/enum.R
$(PACKAGE)/R/enum.R: $(PACKAGE)/src/selfisher.cpp
	echo '## Auto generated - do not edit by hand' > $@
	echo ".valid_link <- c(" >> $@
	grep _link.*= $(PACKAGE)/src/selfisher.cpp | sed s/_link//g >> $@
	echo ")" >> $@

	echo ".valid_family <- c(" >> $@
	grep _family.*= $(PACKAGE)/src/selfisher.cpp | sed s/_family//g >> $@
	echo ")" >> $@

	echo ".valid_covstruct <- c(" >> $@
	grep _covstruct.*= $(PACKAGE)/src/selfisher.cpp | sed s/_covstruct//g >> $@
	echo ")" >> $@

	echo ".valid_ppredictcode <- c(" >> $@
	grep _ppredictcode.*= $(PACKAGE)/src/selfisher.cpp | sed s/_ppredictcode//g >> $@
	echo ")" >> $@

doc-update: $(PACKAGE)/R/*.R
	echo "suppressWarnings(roxygen2::roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\")))" | $(R) --slave
	@touch doc-update

## vignettes

## list of vignette inputs:
rnw_vig += 
rmd_vig += bootstrap catch_comparison_subsampling simple troubleshooting

docdir = $(PACKAGE)/inst/doc
vigdir = $(PACKAGE)/vignettes

docpdf :=  $(rnw_vig:%=${docdir}/%.pdf)
dochtml := $(rmd_vig:%=${docdir}/%.html)

## LaTeX/BibTeX must run in the same directory as the file ...
$(vigdir)/%.pdf: $(vigdir)/%.[Rr]nw
	cd $(vigdir); echo "knitr::knit2pdf(basename(\"$<\"))" | $(R) --slave

%.html: %.rmd
	export NOT_CRAN=true; echo "rmarkdown::render(\"$<\")" | $(R) --slave

## files in doc dir => generate files in vig dir using downstream rules,
## then move them
$(docdir)/%: $(vigdir)/%
	mv $< $@

$(vigdir)/model_evaluation.html: $(vigdir)/model_evaluation.rmd texreg

vignette-update: ${docpdf} ${dochtml}

####
namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(PACKAGE)/R/*.R
	echo "suppressWarnings(roxygen2::roxygenize(\"$(PACKAGE)\",roclets = c(\"namespace\")))" | $(R) --slave

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC)
	$(R) CMD build --resave-data=no $(PACKAGE)

install: $(TARBALL)
	$(R) CMD INSTALL --preclean $<
	@touch $@

## To enable quick compile, run from R:
##    library(TMB); precompile(flags="-O0 -g")
quick-install: enum-update $(PACKAGE)/src/selfisher.so
	$(R) CMD INSTALL $(PACKAGE)

$(PACKAGE)/src/selfisher.so: $(PACKAGE)/src/selfisher.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('selfisher.cpp','-O0 -g')" | $(R) --slave

unexport TEXINPUTS
pdf: $(PACKAGE).pdf
$(PACKAGE).pdf: $(PACKAGE)/man/*.Rd
	rm -f $(PACKAGE).pdf
	$(R) CMD Rd2pdf --no-preview $(PACKAGE)

build:
	$(R) CMD build $(PACKAGE)

check: $(TARBALL)
	$(R) CMD check $(TARBALL)

check-cran: $(TARBALL)
	$(R) CMD check --as-cran $(TARBALL)

## *NOT* using 'R --vanilla' : then cannot find testthat, TMB, etc they are installed into R's "system" library

test:
	echo "devtools::test('selfisher')" | $(R) --slave

quick-check: quick-install ex-test

ex-test:
	echo "library(selfisher); example(selfisher)" | $(R) --slave


unlock:
	\rm -rf `Rscript --vanilla -e 'writeLines(.Library)'`/00LOCK-selfisher
#               ------------------------------------------ = R's system library
#	rm -rf ${R_LIBS}/00LOCK-selfisher
##               ^^^^^^^ This only works if R_LIBS contains a single directory and the same that 'R CMD INSTALL' uses..

clean:
	\rm -f install doc-update
