all: build/v703.pdf

# hier Python-Skripte:
build/charakteristik.pdf: plot.py ../matplotlibrc ../header-matplotlib.tex | build
	# so that matplotlib can find the tex header when running
	# LaTeX in the tmp directory
	# and set the matplotlibrc
	TEXINPUTS=$$(pwd)/..: MATPLOTLIBRC=../matplotlibrc python plot.py

# hier weitere Abhängigkeiten für build/v703.pdf deklarieren:
build/v703.pdf: build/charakteristik.pdf
build/v703.pdf: build/ladung.pdf

build/v703.pdf: FORCE | build
	# to find header and bib files in the main directory
	TEXINPUTS=..: \
	BIBINPUTS=..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	v703.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
