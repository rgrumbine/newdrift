#Shared:
SHARE = math.o metric.o drifter.o io_nc.o

#Offline:
OBJ = $(SHARE) main.o 

#Inline
INLINE = $(SHARE) driver.o mainsub.o

include mk.this

all : drifter in_drifter

drifter : $(OBJ)
	$(FC) $(OBJ)  -L $(NETCDF_LIBRARIES)  -l netcdff -o $(@)

in_drifter : $(INLINE)
	$(FC) $(INLINE)  -L $(NETCDF_LIBRARIES)  -l netcdff -o $(@)

driftlib : $(INLINE)
	ar rv driftlib.a $(INLINE)
	#something for include of mods

clean :
	        rm *.o *.mod 

distclean : clean
	        rm drifter

