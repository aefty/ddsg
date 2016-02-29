ALL: sparsegrid_lib

########################################################################################
# Specify options to compile and create libraries
########################################################################################
CC         = mpicxx
COPTFLAGS  = -O -w
OMAKE      = make
RM         = /bin/rm
AR         = ar
AR_FLAGS   = cr
LIB_SUFFIX = a
RANLIB     = ranlib
#MPICH_DIR  = /opt/mpich-ch_p4-gcc-1.2.7
#MPICH_DIR  = /usr/lib/mpich2/lib
#MPICHPATH = /usr/lib/openmpi/lib

########################################################################################
# compile the source files and create the libraries
########################################################################################

LIB_NAME       = libsparsegrid.$(LIB_SUFFIX)
MAKE_OPTIONS_G = CC="$(CC)" COPTFLAGS="$(COPTFLAGS)" AR="$(AR)" AR_FLAGS="$(AR_FLAGS)" RM="$(RM)"
MAKE_OPTIONS         =   $(MAKE_OPTIONS_G) LIBNAME="$(LIB_NAME)"

########################################################################################
#List all of the source directories you want to compile
########################################################################################

LIBSRCDIRS =  DataStructure AdaptiveDataStructure BlockAllocator AdaptiveSparseGrid Post

sparsegrid_lib:
	-@echo "Begin installing Library"
	for i in $(LIBSRCDIRS); do \
		  (cd $$i && $(OMAKE) lib $(MAKE_OPTIONS)) ; \
	done
	-@echo 
	-@echo "copying libfile to lib "
	-@mv $(LIB_NAME) lib/
	-@echo 
	-@echo "Finish installing Library"
	

clean: 
	$(RM) -f */*.o */*.a

