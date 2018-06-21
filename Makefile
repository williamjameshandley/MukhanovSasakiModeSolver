UNAME=$(shell uname)

ifeq ($(UNAME), Linux)
include Makefile_gnu
else ifeq ($(UNAME), Darwin) 
	include Makefile_osx
endif

source_dir = src
external_dir = external
build_dir = build
binary_dir = bin
lib_dir = $(PWD)/lib

cpp_srcs = $(wildcard $(source_dir)/*.cpp)
cpp_objs = $(cpp_srcs:%.cpp=$(build_dir)/%.o)
cpp_deps = $(cpp_srcs:%.cpp=$(build_dir)/%.d)

f90_srcs = $(wildcard $(source_dir)/*.f90)
f90_objs = $(f90_srcs:%.f90=$(build_dir)/%.o) 


objs = $(cpp_objs) $(f90_objs)

inc += -isystem$(external_dir)/Eigen

#FFLAGS += -J$(build_dir)
libname = msms

LDSHARED = $(LD) -shared

#LDFLAGS += -pthread
LDLIBS += -Wl,-Bstatic -l$(libname) -Wl,-Bdynamic

all: main
extra: tags doc

$(lib_dir)/lib$(libname).so: $(objs)
	@mkdir -p $(@D)
	$(LDSHARED) $^ -o $@ $(LDFLAGS)

$(lib_dir)/lib$(libname).a: $(objs)
	@mkdir -p $(@D)
	$(AR) $(ARFLAGS) $@ $^ 

python $(libdir)/$(libname).so:
	python3 setup.py install --user
	python2 setup.py install --user

# create the directory -- explicitly needed for fortran (annoying)
$(shell mkdir -p $(build_dir)/$(source_dir))

# Compiling the main program
main: $(binary_dir)/main
make_licence: $(binary_dir)/make_licence
main_fortran: $(binary_dir)/main_fortran
lib$(libname): $(lib_dir)/lib$(libname).a
lib$(libname): $(lib_dir)/lib$(libname).so

$(binary_dir)/%: $(build_dir)/%.o lib$(libname)
	@mkdir -p $(@D)
	$(LD) $< -o $@ $(LDFLAGS) -L$(lib_dir) $(LDLIBS) 

# Building general object files and dependencies
$(build_dir)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(inc) -MMD -c $< -o $@

$(build_dir)/%.o: %.f90
	$(FC) $(FFLAGS) $(inc) -c $< -o $@

all_srcs = $(shell find src -name '*.[ch]pp')

# Build tags file
tags: $(all_srcs)
	ctags --extra=+f $(all_srcs)

# Generate documentation
doc: Doxyfile $(all_srcs)
	doxygen $<
	 @echo Documentation in html format may be viewed by opening:
	 @echo "	" file://$(PWD)/doc/html/index.html
	 @echo in a browser

.PHONY: clean main make_licence main_fortran lib$(libname) lib$(libname)_shared

clean:
	$(RM) $(cpp_objs) $(f90_objs) $(cpp_deps) main

purge: clean
	$(RM) -r $(build_dir) $(lib_dir)
	$(RM) -r __pycache__ doc/html
	$(RM) doxygen_warning.log tags 

# Include the dependencies
-include $(cpp_deps)
$(build_dir)/main_fortran.o : $(build_dir)/$(source_dir)/fortran_interface.o
