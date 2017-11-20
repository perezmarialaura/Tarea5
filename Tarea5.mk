*.png: M.dat
	python Plots.py

M.dat:compila
	./a.out

compila: RadialVelocities.dat
	gcc -lm CurvaRotacion.c

clean:
	rm -f M.dat *.png a.out
