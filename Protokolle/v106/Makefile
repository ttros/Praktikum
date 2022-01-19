all: build/v106.pdf

build/v106.pdf: FORCE | build
	# to find header and bib files in the main directory
	TEXINPUTS=..: \
	BIBINPUTS=..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	v106.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
