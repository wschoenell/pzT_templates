source = "paper"

pdf:
	pdflatex $(source)
	bibtex $(source)
	pdflatex $(source)
	pdflatex $(source)

ps:
	latex $(source)
	bibtex $(source)
	latex $(source)
	latex $(source)
	dvips -t a4 -Ppdf $(source) -o $(source).pdf.ps

pdf_from_ps:
	latex $(source)
	bibtex $(source)
	latex $(source)
	latex $(source)
	dvips -t a4 -Ppdf $(source) -o $(source).pdf.ps
	ps2pdf $(source).pdf.ps $(source).pdf

clean:
	rm -f $(source).{aux,lof,log,lot,toc,blg,dvi,tns,synctex.gz,html,4ct,4tc,lg,tmp,xref,css,idv,idx}
	rm -f $(source)_draft.{aux,lof,log,lot,toc,blg,dvi,tns,synctex.gz,html,4ct,4tc,lg,tmp,xref,css,idv,idx}
	rm -f .DS_Store



