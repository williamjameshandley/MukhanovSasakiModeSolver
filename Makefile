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

cpp_src = $(wildcard $(source_dir)/*.cpp)
c_src = $(wildcard $(source_dir)/*.c)
objs = $(cpp_src:%.cpp=$(build_dir)/%.o) $(c_src:%.c=$(build_dir)/%.o)
deps = $(cpp_src:%.cpp=$(build_dir)/%.d) $(c_src:%.c=$(build_dir)/%.d)

inc += -isystem$(external_dir)

libname = msms

LDSHARED = $(LD) -shared

#LDFLAGS += -pthread
LDLIBS +=  -l$(libname)

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

$(binary_dir)/%: $(build_dir)/%.o $(lib_dir)/lib$(libname).so
	@mkdir -p $(@D)
	$(LD) $< -o $@ $(LDFLAGS) -L$(lib_dir) $(LDLIBS) 

# Building general object files and dependencies
$(build_dir)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(inc) -MMD -c $< -o $@

$(build_dir)/%.o: %.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -MMD -c $< -o $@

all_srcs = $(shell find src -name '*.[ch]{,pp}')

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
