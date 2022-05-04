all: build/v500.pdf

# hier Python-Skripte:
build/plot_rot.pdf: plot.py ../matplotlibrc ../header-matplotlib.tex | build
	# so that matplotlib can find the tex header when running
	# LaTeX in the tmp directory
	# and set the matplotlibrc
	TEXINPUTS=$$(pwd)/..: MATPLOTLIBRC=../matplotlibrc python plot.py

# hier weitere Abhängigkeiten für build/v500.pdf deklarieren:
build/v500.pdf: build/plot_rot.pdf
build/v500.pdf: build/plot_gelb.pdf
build/v500.pdf: build/plot_gruen.pdf
build/v500.pdf: build/plot_violett_1.pdf
build/v500.pdf: build/plot_violett_2.pdf

build/v500.pdf: FORCE | build
	# to find header and bib files in the main directory
	TEXINPUTS=..: \
	BIBINPUTS=..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	v500.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
