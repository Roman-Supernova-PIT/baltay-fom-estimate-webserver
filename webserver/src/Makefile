install: snflux snvar prism.py hsiao.txt
	cp -p snflux snvar prism.py hsiao.txt ../

snflux: snflux.f
	gfortran snflux.f -o snflux

snvar: snvar.f
	gfortran snvar.f -o snvar

clean:
	rm -f snvar snflux
	rm -f ../snflux ../snvar ../prism.py ../hsiao.txt

