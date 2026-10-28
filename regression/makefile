
# 1. Compiler and Flag Definitions
FC       = gfortran
FFLAGS   = -O3 -march=native -funroll-loops -ffast-math -std=f2018

# 2. Shared Source Dependencies
# These are the baseline module files required by the test programs.
SHARED_SRCS = math.f90 metric.f90 drifter.f90 kd.f90 kdfind.f90 
SHARED_OBJS = $(SHARED_SRCS:.f90=.o)

SHARED_DIR = ../sorc/shared
MODULES = -I. -I$(SHARED_DIR)
vpath %.f90 $(SHARED_DIR)

# 3. Test Programs (Targets)
# Define all the regression test drivers you want to build.
TEST_SRCS = test_math_utilities.f90 \
            test_kd_mod.f90         \
            test_kdfind.f90 \
            test_metric_matrices.f90  \
            test_drifter_kinematics.f90

# Dynamically generate executable names (e.g., test_kd_mod.f90 -> test_kd_mod)
TARGETS   = $(TEST_SRCS:.f90=)

# ====================================================================
# Build Rules
# ====================================================================

# 'all' builds every single executable listed in $(TARGETS)

.PHONY: all run

all: $(TARGETS)

run: $(TARGETS)
	@echo "trying to run targets"
	@for exec in $(TARGETS); do \
		./$$exec; \
	done

# Pattern Rule for Linking Executables
# This rule tells Make: "To create an executable named X, look for X.o
# and combine it with all the shared object binaries."
%: %.o $(SHARED_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	$@	

# Pattern Rule for Compiling Object Files (.f90 -> .o)
# Applies universally to both shared files and individual test files.
%.o: %.f90
	$(FC) $(FFLAGS) $(MODULES) -c $< -o $@

# ====================================================================
# Explicit Module Dependencies (Crucial for Fortran Compilers)
# ====================================================================
# Fortran requires modules (.mod files) to exist before files that "USE"
# them are compiled. We declare those structural dependencies here:

metric.o: math.o
drifter.o: metric.o math.o
kdfind.o: kd.o
#io_nc.o: drifter.o metric.o

# Map your test drivers to their respective required modules
test_math_utilities.o: math.o
test_kd_mod.o: kd.o
test_kdfind.o: kd.o kdfind.o
test_metric_matrices.o: metric.o
test_drifter_kinematics.o: drifter.o metric.o


#==============================
clean:
	rm *.o *.mod
distclean: clean
	rm $(TARGETS)
