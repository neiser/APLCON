
#       aplcon.a library:  make lib 

FORT        = gfortran
F77OPTIONS  =  -Wunused  -fno-automatic -fno-backslash -O 

lib:
	$(FORT) *.F $(F77OPTIONS)  -c
	ar r /home/volker/libs/aplcon.a  *.o
	rm *.o *~ 
	@echo '... new version of library ~/libs/aplcon.a ready!'



