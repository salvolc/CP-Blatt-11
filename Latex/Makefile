all: main.pdf
	
main.pdf: content/*.tex header.tex lit.bib main.tex | build
	cp lit.bib build/
	cp programme.bib build/
	max_print_line=1048576 TEXINPUTS=build:.: lualatex --output-directory=build --interaction=nonstopmode --halt-on-error main.tex
	BIBINPUTS=build:. biber build/main.bcf
	max_print_line=1048576 TEXINPUTS=build:.: lualatex --output-directory=build --interaction=nonstopmode --halt-on-error main.tex
	max_print_line=1048576 TEXINPUTS=build:.: lualatex --output-directory=build --interaction=nonstopmode --halt-on-error main.tex

build:
	mkdir -p build

clean:
	rm -rf build
pvc: content/*.tex header.tex lit.bib main.tex | build
	cp lit.bib build/
	cp programme.bib build/
	latexmk -lualatex -e '$$pdf_previewer = "evince %O %S";' -interaction=nonstopmode -halt-on-error -synctex=1 --jobname=main -output-directory=build -pvc main.tex