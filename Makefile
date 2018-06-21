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

objs = $(cpp_objs)

inc += -isystem$(external_dir)/Eigen

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

# Compiling the main program
main: $(binary_dir)/main
lib$(libname): $(lib_dir)/lib$(libname).a
lib$(libname): $(lib_dir)/lib$(libname).so

$(binary_dir)/%: $(build_dir)/%.o lib$(libname)
	@mkdir -p $(@D)
	$(LD) $< -o $@ $(LDFLAGS) -L$(lib_dir) $(LDLIBS) 

# Building general object files and dependencies
$(build_dir)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(inc) -MMD -c $< -o $@

all_srcs = $(shell find src -name '*.[ch]pp')

# Build tags file for vim
tags: $(all_srcs)
	ctags --extra=+f $(all_srcs)


.PHONY: clean main lib$(libname)

clean:
	$(RM) $(cpp_objs) $(cpp_deps) main

purge: clean
	$(RM) -r $(build_dir) $(lib_dir)
	$(RM) -r __pycache__ doc/html
	$(RM) tags 

# Include the dependencies
-include $(cpp_deps)
