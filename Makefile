# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/renatoseb/.local/lib/python3.10/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/renatoseb/.local/lib/python3.10/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/renatoseb/2023-2/eda/labs/ss-tree/code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/renatoseb/2023-2/eda/labs/ss-tree/code

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "No interactive CMake dialog available..."
	/home/renatoseb/.local/lib/python3.10/site-packages/cmake/data/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake to regenerate build system..."
	/home/renatoseb/.local/lib/python3.10/site-packages/cmake/data/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/renatoseb/2023-2/eda/labs/ss-tree/code/CMakeFiles /home/renatoseb/2023-2/eda/labs/ss-tree/code//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/renatoseb/2023-2/eda/labs/ss-tree/code/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named ss_tree_test

# Build rule for target.
ss_tree_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 ss_tree_test
.PHONY : ss_tree_test

# fast build rule for target.
ss_tree_test/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_test.dir/build.make CMakeFiles/ss_tree_test.dir/build
.PHONY : ss_tree_test/fast

#=============================================================================
# Target rules for targets named ss_tree_indexing

# Build rule for target.
ss_tree_indexing: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 ss_tree_indexing
.PHONY : ss_tree_indexing

# fast build rule for target.
ss_tree_indexing/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_indexing.dir/build.make CMakeFiles/ss_tree_indexing.dir/build
.PHONY : ss_tree_indexing/fast

#=============================================================================
# Target rules for targets named ss_tree_interface

# Build rule for target.
ss_tree_interface: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 ss_tree_interface
.PHONY : ss_tree_interface

# fast build rule for target.
ss_tree_interface/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/build
.PHONY : ss_tree_interface/fast

#=============================================================================
# Target rules for targets named compile_test

# Build rule for target.
compile_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 compile_test
.PHONY : compile_test

# fast build rule for target.
compile_test/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/compile_test.dir/build.make CMakeFiles/compile_test.dir/build
.PHONY : compile_test/fast

#=============================================================================
# Target rules for targets named compile_indexing

# Build rule for target.
compile_indexing: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 compile_indexing
.PHONY : compile_indexing

# fast build rule for target.
compile_indexing/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/compile_indexing.dir/build.make CMakeFiles/compile_indexing.dir/build
.PHONY : compile_indexing/fast

#=============================================================================
# Target rules for targets named run_test

# Build rule for target.
run_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 run_test
.PHONY : run_test

# fast build rule for target.
run_test/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/run_test.dir/build.make CMakeFiles/run_test.dir/build
.PHONY : run_test/fast

#=============================================================================
# Target rules for targets named run_indexing

# Build rule for target.
run_indexing: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 run_indexing
.PHONY : run_indexing

# fast build rule for target.
run_indexing/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/run_indexing.dir/build.make CMakeFiles/run_indexing.dir/build
.PHONY : run_indexing/fast

#=============================================================================
# Target rules for targets named test

# Build rule for target.
test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 test
.PHONY : test

# fast build rule for target.
test/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test.dir/build.make CMakeFiles/test.dir/build
.PHONY : test/fast

#=============================================================================
# Target rules for targets named indexing

# Build rule for target.
indexing: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 indexing
.PHONY : indexing

# fast build rule for target.
indexing/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/indexing.dir/build.make CMakeFiles/indexing.dir/build
.PHONY : indexing/fast

#=============================================================================
# Target rules for targets named compile_interface

# Build rule for target.
compile_interface: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 compile_interface
.PHONY : compile_interface

# fast build rule for target.
compile_interface/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/compile_interface.dir/build.make CMakeFiles/compile_interface.dir/build
.PHONY : compile_interface/fast

#=============================================================================
# Target rules for targets named run_interface

# Build rule for target.
run_interface: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 run_interface
.PHONY : run_interface

# fast build rule for target.
run_interface/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/run_interface.dir/build.make CMakeFiles/run_interface.dir/build
.PHONY : run_interface/fast

#=============================================================================
# Target rules for targets named interface

# Build rule for target.
interface: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 interface
.PHONY : interface

# fast build rule for target.
interface/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/interface.dir/build.make CMakeFiles/interface.dir/build
.PHONY : interface/fast

CortexAPI.o: CortexAPI.cpp.o
.PHONY : CortexAPI.o

# target to build an object file
CortexAPI.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/CortexAPI.cpp.o
.PHONY : CortexAPI.cpp.o

CortexAPI.i: CortexAPI.cpp.i
.PHONY : CortexAPI.i

# target to preprocess a source file
CortexAPI.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/CortexAPI.cpp.i
.PHONY : CortexAPI.cpp.i

CortexAPI.s: CortexAPI.cpp.s
.PHONY : CortexAPI.s

# target to generate assembly for a file
CortexAPI.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/CortexAPI.cpp.s
.PHONY : CortexAPI.cpp.s

Interface.o: Interface.cpp.o
.PHONY : Interface.o

# target to build an object file
Interface.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/Interface.cpp.o
.PHONY : Interface.cpp.o

Interface.i: Interface.cpp.i
.PHONY : Interface.i

# target to preprocess a source file
Interface.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/Interface.cpp.i
.PHONY : Interface.cpp.i

Interface.s: Interface.cpp.s
.PHONY : Interface.s

# target to generate assembly for a file
Interface.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/Interface.cpp.s
.PHONY : Interface.cpp.s

SStree.o: SStree.cpp.o
.PHONY : SStree.o

# target to build an object file
SStree.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_test.dir/build.make CMakeFiles/ss_tree_test.dir/SStree.cpp.o
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_indexing.dir/build.make CMakeFiles/ss_tree_indexing.dir/SStree.cpp.o
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/SStree.cpp.o
.PHONY : SStree.cpp.o

SStree.i: SStree.cpp.i
.PHONY : SStree.i

# target to preprocess a source file
SStree.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_test.dir/build.make CMakeFiles/ss_tree_test.dir/SStree.cpp.i
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_indexing.dir/build.make CMakeFiles/ss_tree_indexing.dir/SStree.cpp.i
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/SStree.cpp.i
.PHONY : SStree.cpp.i

SStree.s: SStree.cpp.s
.PHONY : SStree.s

# target to generate assembly for a file
SStree.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_test.dir/build.make CMakeFiles/ss_tree_test.dir/SStree.cpp.s
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_indexing.dir/build.make CMakeFiles/ss_tree_indexing.dir/SStree.cpp.s
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/SStree.cpp.s
.PHONY : SStree.cpp.s

indexing.o: indexing.cpp.o
.PHONY : indexing.o

# target to build an object file
indexing.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_indexing.dir/build.make CMakeFiles/ss_tree_indexing.dir/indexing.cpp.o
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/indexing.cpp.o
.PHONY : indexing.cpp.o

indexing.i: indexing.cpp.i
.PHONY : indexing.i

# target to preprocess a source file
indexing.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_indexing.dir/build.make CMakeFiles/ss_tree_indexing.dir/indexing.cpp.i
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/indexing.cpp.i
.PHONY : indexing.cpp.i

indexing.s: indexing.cpp.s
.PHONY : indexing.s

# target to generate assembly for a file
indexing.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_indexing.dir/build.make CMakeFiles/ss_tree_indexing.dir/indexing.cpp.s
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/indexing.cpp.s
.PHONY : indexing.cpp.s

main.o: main.cpp.o
.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_test.dir/build.make CMakeFiles/ss_tree_test.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i
.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_test.dir/build.make CMakeFiles/ss_tree_test.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s
.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_test.dir/build.make CMakeFiles/ss_tree_test.dir/main.cpp.s
.PHONY : main.cpp.s

tinyfiledialogs.o: tinyfiledialogs.c.o
.PHONY : tinyfiledialogs.o

# target to build an object file
tinyfiledialogs.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/tinyfiledialogs.c.o
.PHONY : tinyfiledialogs.c.o

tinyfiledialogs.i: tinyfiledialogs.c.i
.PHONY : tinyfiledialogs.i

# target to preprocess a source file
tinyfiledialogs.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/tinyfiledialogs.c.i
.PHONY : tinyfiledialogs.c.i

tinyfiledialogs.s: tinyfiledialogs.c.s
.PHONY : tinyfiledialogs.s

# target to generate assembly for a file
tinyfiledialogs.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ss_tree_interface.dir/build.make CMakeFiles/ss_tree_interface.dir/tinyfiledialogs.c.s
.PHONY : tinyfiledialogs.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... compile_indexing"
	@echo "... compile_interface"
	@echo "... compile_test"
	@echo "... indexing"
	@echo "... interface"
	@echo "... run_indexing"
	@echo "... run_interface"
	@echo "... run_test"
	@echo "... test"
	@echo "... ss_tree_indexing"
	@echo "... ss_tree_interface"
	@echo "... ss_tree_test"
	@echo "... CortexAPI.o"
	@echo "... CortexAPI.i"
	@echo "... CortexAPI.s"
	@echo "... Interface.o"
	@echo "... Interface.i"
	@echo "... Interface.s"
	@echo "... SStree.o"
	@echo "... SStree.i"
	@echo "... SStree.s"
	@echo "... indexing.o"
	@echo "... indexing.i"
	@echo "... indexing.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... tinyfiledialogs.o"
	@echo "... tinyfiledialogs.i"
	@echo "... tinyfiledialogs.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

