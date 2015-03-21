
#

FORT       = gfortran
F77OPTIONS = -Wunused  -fno-automatic -fno-backslash -O 
RLIB       = /lib/libreadline.so.6

# create new load-module 
#                  under linux

avmain:   avmain.F avlass.F avvsys.F  
	$(FORT) avmain.F  \
	avlass.F avvsys.F \
	$(F77OPTIONS) -o avmain  \
	/home/volker/libs/aplcon.a \
	$(RLIB) -lncurses
	@echo '... new version of avmain ready - try it out!'


########   clean

.PHONY: clean
clean:
	-rm  avmain  *~  *.o 

