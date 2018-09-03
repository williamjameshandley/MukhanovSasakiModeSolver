UNAME=$(shell uname)

ifeq ($(UNAME), Linux)
include Makefile_gnu
else ifeq ($(UNAME), Darwin) 
	include Makefile_osx
endif

source_dir = $(PWD)/src
test_dir = $(source_dir)/test
external_dir = $(PWD)/external
cephes_dir = $(external_dir)/cephes
odepack_dir = $(external_dir)/odepack
build_dir = $(PWD)/build
binary_dir = $(PWD)/bin
lib_dir = $(PWD)/lib

makefiles = $(wildcard Makefile*)

src = $(wildcard $(source_dir)/*.cpp)
cephes_src := $(wildcard $(cephes_dir)/*.c) 
odepack_cpp_src := $(wildcard $(odepack_dir)/*.cpp)
odepack_f_src := $(wildcard $(odepack_dir)/*.f)

test_srcs = $(wildcard $(test_dir)/*.cpp)
test_objs = $(test_srcs:%.cpp=$(build_dir)/%.o)
test_deps = $(test_srcs:%.cpp=$(build_dir)/%.d)

objs = $(src:%.cpp=$(build_dir)/%.o) $(cephes_src:%.c=$(build_dir)/%.o) $(odepack_f_src:%.f=$(build_dir)/%.o) $(odepack_cpp_src:%.cpp=$(build_dir)/%.o) 
deps = $(src:%.cpp=$(build_dir)/%.d) $(cephes_src:%.c=$(build_dir)/%.d)

inc += -isystem$(external_dir) -I$(PWD)/src

libname = msms

LDSHARED = $(LD) -shared

#LDFLAGS += -pthread
LDLIBS +=  -l$(libname)

all: test main
extra: tags doc
main: $(binary_dir)/main 
test: $(binary_dir)/test 

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
$(binary_dir)/main: $(build_dir)/main.o $(objs)
	@mkdir -p $(@D)
	$(LD) $^ -o $@ $(LDFLAGS)

$(binary_dir)/test: $(build_dir)/test_main.o $(objs) $(test_objs)
	@mkdir -p $(@D)
	$(LD) $^ -o $@ $(LDFLAGS)

# Building general object files and dependencies
$(build_dir)/%.o: %.cpp $(makefiles)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(inc) -MMD -c $< -o $@

$(build_dir)/%.o: %.c $(makefiles)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(inc) -MMD -c $< -o $@

$(build_dir)/%.o: %.f $(makefiles)
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) -c $< -o $@

all_srcs = $(shell find src -name '*.[ch]pp')

# Build tags file for vim
tags: $(all_srcs)
	ctags --extra=+f $(all_srcs)


.PHONY: clean main

clean:
	$(RM) $(objs) $(deps) main

purge: clean
	$(RM) -r $(build_dir) $(lib_dir)
	$(RM) -r __pycache__ doc/html
	$(RM) tags 

# Include the dependencies
-include $(deps)
$(build_dir)/external/odepack/opkda1.o:
$(build_dir)/external/odepack/opkda2.o:
$(build_dir)/external/odepack/opkdmain.o:
