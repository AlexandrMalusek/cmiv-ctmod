CTMOD = ../bin/ctmod

output.txt : input.txt
	$(CTMOD) < $< > $@

symlinks :
	ln -fs $(CTMODSYS)/src/obj obj
	ln -fs $(CTMODSYS)/material material
	ln -fs $(CTMODSYS)/detector detector 
	ln -fs $(CTMODSYS)/spectra spectra
