#Shared:
SHARE = math.o metric.o drifter.o io_nc.o

#Offline:
OBJ = $(SHARE) main.o 

#Inline
INLINE = $(SHARE) read.o write.o run.o 

include mk.this

drifter : $(OBJ)
	echo zzz netcdf_lib = $(NETCDF_LIBRARIES)
	$(FC) $(OBJ)  -L $(NETCDF_LIBRARIES)  -l netcdff -o $(@)

driftlib : $(INLINE)
	ar rv driftlib.a $(INLINE)
	#something for include of mods

clean :
	        rm *.o *.mod 

distclean : clean
	        rm drifter

