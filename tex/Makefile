# Make file for creating .pdf file from .tex file
# ** always "make clean" first!
# Originator: RBM, Oct'05
#
# detect the machine architecture: Linux/Mac OS X 
machine := ${shell uname -m}
#
# Mac OS X
#ifeq (${machine}, i386)
LATEX=pdflatex
BIBTEX=bibtex
DVIPDF=dvipdf 
#ACROBAT=acroread
ACROBAT=kpdf

#endif
# Linux 
#ifeq (${machine}, i686)
#LATEX=/usr/bin/pdflatex
#BIBTEX=/usr/bin/bibtex
#ACROBAT=/usr/bin/acroread
#endif
#
#------------------------------------------------------------
# remove ("clean") object files:
RM = rm -f
#
#------------------------------------------------------------
# define source, object and target files
NAME = thesis
DIR = tex
SOURCES = $(NAME).tex 
OBJECTS = *.aux *.bbl *.blg *.lof *.lot *.log *.toc
TARGETS = $(NAME).pdf
DVI = $(NAME).dvi
#
#---------------------- make --------------------------------
#all : clean thesis pdf clean
all : clean thesis clean
#
compile : 
	$(LATEX) $(SOURCES)
#
dvi : clean thesis clean
#
thesis : 
	$(LATEX) $(SOURCES); \
	$(BIBTEX) $(NAME); \
	$(BIBTEX) $(NAME); \
	$(LATEX) $(SOURCES); \
	$(LATEX) $(SOURCES)   
#
#pdf :
#	$(DVIPDF) $(DVI); \
#	$(RM) $(OBJECTS)   
#	$(ACROBAT) $(TARGETS) & 
#
clean:
	$(RM) $(OBJECTS)
#------------------------------------------------------------
