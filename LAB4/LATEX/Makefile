all: report.pdf

report.pdf: main.tex
	texi2pdf main.tex -o report.pdf

run: report.pdf
	killall zathura    || true
	zathura report.pdf &

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb

ultraclean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb report.pdf

