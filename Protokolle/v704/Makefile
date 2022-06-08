all: build/v704.pdf

# hier Python-Skripte:
build/beta.pdf: plot.py ../matplotlibrc ../header-matplotlib.tex | build
	# so that matplotlib can find the tex header when running
	# LaTeX in the tmp directory
	# and set the matplotlibrc
	TEXINPUTS=$$(pwd)/..: MATPLOTLIBRC=../matplotlibrc python plot.py

# hier weitere Abhängigkeiten für build/v704.pdf deklarieren:
build/v704.pdf: build/beta.pdf
# build/v704.pdf: build/Rhodium.pdf

build/v704.pdf: FORCE | build
	# to find header and bib files in the main directory
	TEXINPUTS=..: \
	BIBINPUTS=..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	v704.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
