BUILD  := build
ROOTDIR := .

# This line will check if HEADAS is set, and otherwise will fall back to using
# Python and xspectrampoline to configure it:
HEADAS ?= $(shell python3 -c "import xspectrampoline_helpers as h ; print(h.get_HEADAS())")
HEADAS_LIB := ${HEADAS}/lib
HEADAS_INCLUDE := ${HEADAS}/include

VERSION := $(shell cat VERSION)

# A cache directory for static files that should be persistent between `make
# clean` calls.
CACHEDIR := cache
# Where the reltrans tables are stored. This will use the user's environment
# variable if set.
RELTRANS_TABLES ?= $(CACHEDIR)/tables

# These may be set when invoking `make`, such as `make DEBUG=1 SANITIZE=1`.
# The `DEBUG` option compiles a debug build of reltrans (see below).
DEBUG = 0
# The `SANITIZE` option enables the address sanitizer in the library (see
# below).
SANITIZE = 0
# The `COVERAGE` option compiles the library with coverage enabled. This can be
# used with the `cov-report` target to generate a coverage breakdown.
COVERAGE = 0

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
		  -I$(BUILD)/include \
		  -I$(HEADAS_INCLUDE) \
		  -I$(HEADAS_INCLUDE)/fftw \
		  -J$(BUILD)/cache \
		  -I$(BUILD)/cache

LDFLAGS := -lkerrz -Wl,-rpath,'$(abspath $(BUILD)/lib)' -L$(BUILD)/lib \
	-L$(HEADAS_LIB) -lXSFunctions -lXSModel -lfftw3 \
	$(shell ls -1 $(HEADAS_LIB)/libcfitsio.*$(SHARED_EXT)* | head -n1) \
	-Wl,-rpath,'$(HEADAS_LIB)'

# Only pass the version macro if git found a tag
ifneq ($(VERSION),)
	FFLAGS += -DRELTRANS_VERSION='"$(VERSION)"'
endif

ifeq ($(DEBUG),1)
	# Compile reltrans in 'debug' mode, which means disabling optimisations and
	# including debug symbols.
	FFLAGS += -DRELTRANS_BUILD_DEBUG=1 -g
	CFLAGS += -g
else
	# The default arguments used to compile reltrans
	FFLAGS += -O3
	CFLAGS += -O3
	# For production releases also enable link-time optimisations
	LDFLAGS += -flto
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

ifeq ($(COVERAGE),1)
	# `-coverage` should be expanded to the other two but it isn't for all
	# Fortran compiler versions.
	FFLAGS += -coverage -fprofile-arcs -ftest-coverage
endif

ifeq ($(TARGET),Linux)
	FFLAGS += -shared -export-dynamic
	LDFLAGS += -lm -lpthread
	LDRPATH := -Wl,-rpath,'$$ORIGIN/../lib'
	# These flags are only passed to the executables that are not directly part
	# of the reltrans library.
	EXE_LDFLAGS := -lmvec
	SHARED_EXT := so
	SED_INPLACE = sed -i
else
ifeq ($(TARGET),Darwin)
	FFLAGS += -dynamiclib
	LDFLAGS += -lgfortran
	LDRPATH := -Wl,-rpath,@executable_path/../lib
	EXE_LDFLAGS :=
	SHARED_EXT := dylib
	# MacOS sed needs an extra useless argument
	SED_INPLACE = sed -i ''
endif
endif

# the path to the reltrans library for the -L linker flag
LIB_PATH := $(abspath $(BUILD)/lib)
RELTRANS_SHARED_LIBRARY := $(BUILD)/lib/libreltrans.$(SHARED_EXT)
LIB_KERRZ = $(shell python3 -c "import kerrz_lib ; print(kerrz_lib.bindings.KERRZ_PATH)")
LIB_KERRZ_SYMLINK := $(BUILD)/lib/libkerrz.$(SHARED_EXT)
KERRZ_F90 := $(BUILD)/cache/kerrz.f90

all: $(BUILD) $(RELTRANS_SHARED_LIBRARY)

exe: $(BUILD)/bin/relcli

dummy: $(BUILD)/bin/dummy

.PHONY: coverage
coverage: coverage.info

coverage.info:
	# Requires both `gcov` and `lcov`
	gcov -o $(BUILD)/cache -b wrappers.f90
	lcov --gcov-tool gcov --capture --directory . --output-file coverage.info

.PHONY: covreport
covreport: coverage.info
	genhtml --output-directory $(BUILD)/html coverage.info

$(BUILD)/bin/relcli: ./utils/cli.c $(BUILD)/lib/libreltrans.$(SHARED_EXT)
	$(CC) $(CFLAGS) utils/cli.c -o $@ \
		-L$(BUILD)/lib -lgfortran -lc -lm $(EXE_LDFLAGS) \
		$(LDRPATH) -lreltrans

$(BUILD)/bin/dummy: ./utils/dummy.c $(BUILD)/lib/libreltrans.$(SHARED_EXT)
	$(CC) $(CFLAGS) utils/dummy.c -o $@ \
		-L$(BUILD)/lib -lgfortran -lc -lm $(EXE_LDFLAGS) \
		$(LDRPATH) -lreltrans

# Need to use abspath here so that on MacOS the correct linker identity is
# generated. Macos does library pathing differently, and the easiest thing to
# do is to make sure anything that links against reltrans gets the absolute
# path to the library.
#
# Note: this does mean that the library cannot be relocated on the machine. If
# someone wants to install it to a different location, the easiest thing to do
# would be to either tell them to run `make BUILD=/path/to/opt/`, or to invoke
# `install_name_tool` (see discussion in PR #55).
$(RELTRANS_SHARED_LIBRARY): $(BUILD)/cache/wrappers.o $(BUILD)/cache/kerrz.o $(BUILD)/cache/constants.o
	$(FC) $(FFLAGS) $^ -o $(abspath $@) $(LDFLAGS)

$(BUILD)/cache/wrappers.o: $(ROOTDIR)/wrappers.f90 $(BUILD)/cache/kerrz.o $(BUILD)/cache/constants.o $(ALL_RELTRANS_SOURCE_FILES)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/cache/%.o: $(ROOTDIR)/subroutines/%.f90
	mkdir -p $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)/bin $(BUILD)/lib $(BUILD)/include $(BUILD)/share $(BUILD)/cache

.PHONY: clean
clean:
	rm -rf $(BUILD)
	rm -rf *.f90.gcov
	rm -rf coverage.info

.PHONY: format
format:
	clang-format -i ./utils/cli.c

$(BUILD)/cache/kerrz.o: $(KERRZ_F90) $(LIB_KERRZ_SYMLINK)
	$(FC) $(FFLAGS) -c $< -o $@

$(LIB_KERRZ_SYMLINK):
	# Assume kerrz is pip-installed
	# so just symlink the library so we have the library in a convenient place
	ln -s "$(LIB_KERRZ)" "$@"

$(KERRZ_F90):
	ln -s $(shell python3 -c "import kerrz_lib ; print(kerrz_lib.bindings.KERRZ_PATH.parent / 'kerrz.f90')") "$@"

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

.PHONY: tables-renorm
tables-renorm:
	# Normalise the tables
	python3 ./renormalise_table.py

.PHONY: tables-fetch
fetch-tables: $(RELTRANS_TABLES)
	@echo "Tables located at '$(RELTRANS_TABLES)'"
	@echo "Please run"
	@echo ""
	@echo "    export RELTRANS_TABLES=$(RELTRANS_TABLES)"
	@echo ""
	@echo "to instruct reltrans to use that path. Add to your `~/.bashrc` to"
	@echo "make the change persistent. To redownload the tables, use"
	@echo ""
	@echo "    unset RELTRANS_TABLES"
	@echo ""
	@echo "and remove the '$(CACHEDIR)/tables' directory."

$(CACHEDIR)/tables:
	mkdir -p $@
	@echo "Downloading pre-normalised tables (may take a few minutes)..."
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillver-a-Ec5_normalised.fits \
		-o $(@)/xillver-a-Ec5_normalised.fits
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillverCp_v3.4_normalised.fits-00 \
		-o $(@)/xillverCp_v3.4_normalised.fits-00
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillverCp_v3.4_normalised.fits-01 \
		-o $(@)/xillverCp_v3.4_normalised.fits-01
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillverCp_v3.4_normalised.fits-02 \
		-o $(@)/xillverCp_v3.4_normalised.fits-02
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillverCp_v3.4_normalised.fits-03 \
		-o $(@)/xillverCp_v3.4_normalised.fits-03
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillverCp_v3.4_normalised.fits-04 \
		-o $(@)/xillverCp_v3.4_normalised.fits-04
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillverCp_v3.4_normalised.fits-05 \
		-o $(@)/xillverCp_v3.4_normalised.fits-05
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillverCp_v3.4_normalised.fits-06 \
		-o $(@)/xillverCp_v3.4_normalised.fits-06
	curl -sL \
		https://github.com/reltrans/model-data/releases/download/v0.1.0/xillverD-5_normalised.fits \
		-o $(@)/xillverD-5_normalised.fits
	( cd $(CACHEDIR)/tables && \
		cat `ls xillverCp_v3.4_normalised.fits-* | sort -V` > xillverCp_v3.4_normalised.fits )

.PHONY: instrument-files
instrument-files:
	@echo "Downloading test suite instrument files..."
	mkdir -p $(CACHEDIR)/instrument-files
	curl -sL \
		"https://github.com/reltrans/model-data/releases/download/v0.1.0/nicer-consim135p-teamonly-array50.arf" \
		-o $(CACHEDIR)/instrument-files/nicer-consim135p-teamonly-array50.arf
	curl -sL \
		"https://github.com/reltrans/model-data/releases/download/v0.1.0/nicer-rmf6s-teamonly-array50.rmf" \
		-o $(CACHEDIR)/instrument-files/nicer-rmf6s-teamonly-array50.rmf
	@echo "Instrument files now located at '$(CACHEDIR)/instrument-files'"
