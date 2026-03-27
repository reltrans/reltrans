BUILD  := build
ROOTDIR := .
HEADAS_LIB := ${HEADAS}/lib
HEADAS_INCLUDE := ${HEADAS}/include

# Derive the version from the most recent git tag, or from .git_archival.txt
# (automatically populated by git archive / GitHub release tarballs).
VERSION := $(shell git fetch --tags --quiet 2>/dev/null; git describe --tags --always 2>/dev/null || sed -n 's/describe-name: *//p' .git_archival.txt 2>/dev/null | grep -v 'Format:')
VERSION := $(shell echo '$(VERSION)' | sed 's/-.*//')


# These may be set when invoking `make`, such as `make DEBUG=1 SANITIZE=1`.
# The `DEBUG` option compiles a debug build of reltrans (see below).
DEBUG = 0
# The `SANITIZE` option enables the address sanitizer in the library (see
# below).
SANITIZE = 0

# The name to use for the reltrans library (must be different from reltrans, as
# libreltrans is the compiled reltrans library)
XSPEC_RELTRANS_NAME = xsreltrans

# This are configurable from the command line
TARGET = $(shell uname)
FC = gfortran

ALL_RELTRANS_SOURCE_FILES = $(shell find $(ROOTDIR)/subroutines -name '*.f90')

CFLAGS := -fno-omit-frame-pointer

FFLAGS := -cpp -DHAVE_INLINE \
		  -fPIC -fno-automatic -fno-second-underscore \
		  -fno-omit-frame-pointer \
		  -fopenmp \
		  -I$(BUILD)/include \
		  -I$(HEADAS_INCLUDE) \
		  -I$(HEADAS_INCLUDE)/fftw \
		  -J$(BUILD)/cache \
		  -I$(BUILD)/cache

# Only pass the version macro if git found a tag
ifneq ($(VERSION),)
FFLAGS += -DRELTRANS_VERSION='"$(VERSION)"'
endif

LDFLAGS := -lXSFunctions -lXSModel -lfftw3 -lcfitsio \
		   -L$(BUILD)/lib -L$(HEADAS_LIB)

ifeq ($(DEBUG),1)
	# Compile reltrans in 'debug' mode, which means disabling optimisations and
	# including debug symbols.
	FFLAGS += -g
	CFLAGS += -g
else
	# The default arguments used to compile reltrans
	FFLAGS += -O3
	CFLAGS += -O3
endif

ifeq ($(SANITIZE),1)
	# Include the address sanitizer. This is a compiler feature which adds a
	# runtime address sanitizer that checks whether all memory addresses being
	# accessed are valid (i.e. avoiding buffer overflows, use-after-frees, and
	# so on).
	#
	# This has a runtime overhead, so must be explicitly enabled. If the
	# sanitizer is triggered, it will print a traceback and some memory
	# information, making it easier to locate the line which accessed illegal
	# memory.
	FFLAGS += -fsanitize=address
	CFLAGS += -fsanitize=address
endif

ifeq ($(TARGET),Linux)
	FFLAGS += -shared -export-dynamic
	LDFLAGS += -lm -lpthread
	SHARED_EXT := so
	SED_INPLACE = sed -i
else
ifeq ($(TARGET),Darwin)
	FFLAGS += -dynamiclib
	LDFLAGS += -lgfortran
	SHARED_EXT := dylib
	# MacOS sed needs an extra useless argument
	SED_INPLACE = sed -i ''
endif
endif

# the path to the reltrans library for the -L linker flag
LIB_PATH := $(abspath $(BUILD)/lib)
RELTRANS_SHARED_LIBRARY := $(BUILD)/lib/libreltrans.$(SHARED_EXT)

all: $(BUILD) $(RELTRANS_SHARED_LIBRARY)

exe: $(BUILD)/bin/relcli

$(BUILD)/bin/relcli: ./utils/cli.c $(BUILD)/lib/libreltrans.$(SHARED_EXT)
	$(CC) $(CFLAGS) utils/cli.c -o $(BUILD)/bin/relcli \
		-L$(BUILD)/lib -lgfortran -lc -lm -lmvec \
		-Wl,-rpath,'$$ORIGIN/../lib' -lreltrans

# Need to use abspath here so that on MacOS the correct linker identity is
# generated. Macos does library pathing differently, and the easiest thing to
# do is to make sure anything that links against reltrans gets the absolute
# path to the library.
#
# Note: this does mean that the library cannot be relocated on the machine. If
# someone wants to install it to a different location, the easiest thing to do
# would be to either tell them to run `make BUILD=/path/to/opt/`, or to invoke
# `install_name_tool` (see discussion in PR #55).
$(RELTRANS_SHARED_LIBRARY): $(BUILD) $(BUILD)/cache/wrappers.o
	$(FC) $(FFLAGS) $(BUILD)/cache/wrappers.o -o $(abspath $@) $(LDFLAGS)

$(BUILD)/cache/%.o: $(ROOTDIR)/%.f90 $(ALL_RELTRANS_SOURCE_FILES)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)/bin $(BUILD)/lib $(BUILD)/include $(BUILD)/share $(BUILD)/cache

.PHONY: clean
clean:
	rm -rf $(BUILD)

.PHONY: format
format:
	clang-format -i ./utils/cli.c

.PHONY: xspec
xspec: $(RELTRANS_SHARED_LIBRARY) xspec/lmodel_reltrans.dat xspec/compile_reltrans.xcm
	# Copy the necessary XSPEC files into the build directory
	cp -r xspec $(BUILD)/xspec
	cd $(BUILD)/xspec && \
		echo "initpackage $(XSPEC_RELTRANS_NAME) lmodel_reltrans.dat .\n exit" | xspec
	# Delete the immediately compiled shared library, so we can recompile it
	# with our options.
	rm $(BUILD)/xspec/libxsreltrans.$(SHARED_EXT)
	# Patch the XSPEC generated Makefile so that it uses the shared library
	# compiled outside of XSPEC
	$(SED_INPLACE) 's|-lXSFunctions|-lXSFunctions -L$(LIB_PATH) -Wl,-rpath,"$(LIB_PATH)" -lreltrans|g' \
		$(BUILD)/xspec/Makefile
	# Set the library name
	$(SED_INPLACE) 's|{LIBRARY_NAME}|$(XSPEC_RELTRANS_NAME)|g' $(BUILD)/xspec/compile_reltrans.xcm
	# Compile and pray XSPEC is happy
	cd $(BUILD)/xspec && xspec - compile_reltrans.xcm
	@echo "--------------------------------------------------------------------"
	@echo "Build succeeded and all XSPEC checks passed."
	@echo ""
	@echo "To use reltrans in XSPEC, start XSPEC and load the model:"
	@echo ""
	@echo "    lmod xsreltrans $(abspath $(BUILD)/xspec/)"
	@echo ""
	@echo "For more information, consult the reltrans documentation (see the"
	@echo "README included in the repository)."

.PHONY: tables
tables:
	# Normalise the tables
	python3 ./renormalise_table.py
