SHELL=/bin/sh

# view path for 'make'
# --------------------
# 'src' is defined in the shell environment
VPATH=../src

# compiler and options
# --------------------
# include fortran compiler command and options
include Makeinclude_fc

# suffix rules
# ------------
include ../src/Makeinclude_suffix

# dependencies not covered by suffix rules
# and list of executables (${programs})
# ----------------------------------------
EXT_SUFFIX := $(shell python3-config --extension-suffix)
include ../src/Makeinclude_deps

# utilities
# ---------
lib : ${lib}
programs: ${lib} ${programs}
wrapper: ${lib} ${wrapper}
all: ${lib} ${programs} ${wrapper}
clean    :
	rm -f *.o *.so *.${MODEXT} ${programs}
	cd ../bin; rm -f ${programs} *.nc *.so
