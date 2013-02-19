SUBDIRS = src test examples

EXTRA_DIST = Doxyfile html latex/refman.pdf

documentation :
	doxygen
	cd latex && pdflatex refman.tex

