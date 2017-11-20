*.pdf:*.png
	pdflatex Results_hw5.tex
*.png:M.dat
	python Plots.py
M.dat:compila
	./a.out
compila:RadialVelocities.dat
	gcc CurvaRotacion.c -lm
clean:
	rm -f M.dat *.png a.out *.aux *.log *.pdf
