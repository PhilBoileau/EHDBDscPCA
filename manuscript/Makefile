MAIN = main
RESP = response

.PHONY : all
all : doc clean

doc :
	pdflatex $(MAIN).tex
	bibtex $(MAIN)
	pdflatex $(MAIN).tex
	pdflatex $(MAIN).tex

resp:
	pdflatex $(RESP).tex
	bibtex $(RESP)
	pdflatex $(RESP).tex
	pdflatex $(RESP).tex

docminimal:
	pdflatex $(MAIN).tex

clean :
	rm -f $(addprefix $(MAIN),.aux .blg .fff .log .mw .out)
	rm -f $(addprefix $(RESP),.aux .blg .fff .log .mw .out)
