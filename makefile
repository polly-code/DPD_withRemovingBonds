objects = ChromDPDr.f90

   FC = mpif90		
   FFLAGS = 
   LDFLAGS = -O2 -ffree-line-length-512
   TARGET = dpd

default: $(objects) 
	$(FC) $(LDFLAGS) -o $(TARGET) $(objects)
   $(objects) :

clean: 
	rm -f $(TARGET) *.o
