#---------------------------C++ Makefile--------------------------------------#
# 
# Features: Handles a multi-directory structure for projects and automatic
# component-specific compilation for dependencies. It also uses an iterated 
# style of recursive make to build components in order. In addition, if
# needed source files are outside of the project directory, there are options
# to draw them in also.
# 
# Type make to compile, make clean to remove executable and object files,
# and make cleaner to remove TARGET, DDIR, LDIR, and ODIR. Using this make file 
# in an empty directory will initialize the folder hierarchy.
# 
# Todo: Make compatible with .c files.
# 
# Layout:
#	.		: Put this Makefile in the top level of the project
#	./bin		: Directory for where to place the target executable.
#	./obj		: Directory for object file storage.
#	./.obj/.dep	: Directory for placing dependency information.
#	./.obj/.lnk	: Directory for placing links to source code.
#	./src		: Directory for source files. Subdirectories of ./src
#			  are automatically drawn in via to search for source files
#                         through use of the find tool.
# 
#    See http://www.gnu.org/software/make/manual/make.html for more info.      #

####################        CHANGEABLE ITEMS     ###############################

#------------------  Output File and Folder Hierarchy -------------------------#

# Change the name of the compiled program here.
NAME = myprogram

# If desired, change these variables to store files in other locations.
TDIR = ./bin
ODIR = ./.obj
LDIR = $(ODIR)/.lnk
DDIR = $(ODIR)/.dep
SDIR = ./src

# Use EXTDIRS to locate external directories with required source files.
# Make sure that the name of every file pulled into the project is unique. 
# If the source file has spaces in the path, replace them with '\ '.
# Note that subdirectories of paths in EXTDIRS are not searched for.

EXTDIRS = # ex. ./myproject/code/ ./other\ /project/code

####################     DON'T CHANGE BELOW HERE     ###########################

#-------------------  Variable Defintions -------------------------------------#

# Flags for the compiler.
CC = g++
CPPFLAGS = -g -Wall -Werror -std=c++17
DEPFLAGS = -MT $@ -MM -MP -MF # when using $(DEPFLAGS), put output info after.

sp :=\\

# SEDIR is Source + External directories
SEDIR =

# Name and location of final executable.
$(shell if [ ! -d "$(TDIR)/" ]; then mkdir $(TDIR)/; fi)
TARGET := $(TDIR)/$(NAME)

# Grab source files from $(SDIR). Look through ./src and directories named in $(EXTDIRS)
INIT := $(shell ! test -d $(SDIR) && echo init)
$(shell if [ ! -d "$(SDIR)/" ]; then mkdir $(SDIR)/; fi)
SEDIR := $(shell find $(SDIR) -type d -print | sed 's, ,\ ,' | sed 's,\(.*\),"\1",')
SEDIR += $(EXTDIRS)
SRC := $(shell find $(SEDIR) -maxdepth 1 -name '*.cpp' -print | sed 's,.*/\(.*\)\.cpp,\1,')

# Convert the source file names to object file names.
$(shell if [ ! -d "$(ODIR)/" ]; then mkdir $(ODIR)/; fi)
OBJS := $(foreach file, $(SRC), $(ODIR)/$(file).o)

# Process source file information to set up links directories.
$(shell if [ ! -d "$(LDIR)/" ]; then mkdir $(LDIR)/; fi)
LNKPATH := $(foreach file, $(SRC), $(LDIR)/$(file)-dir)
vpath %.cpp $(LNKPATH)

# Prepare dependency information.
$(shell if [ ! -d "$(DDIR)/" ]; then mkdir $(DDIR)/; fi)
DEPS := $(foreach file, $(SRC), $(DDIR)/$(file).d)

# Prepare switches to handle the execution of certain parts each time MAKE is called.
STEPS = ZerothStep FirstStep SecondStep ThirdStep NullStep InitStep
TAGS = Tag0 Tag1 Tag2 Tag3 TagNull TagInit

ifeq ($(MAKELEVEL), 0)
STEP = ZerothStep
endif
ifeq ($(MAKELEVEL), 1)
STEP = FirstStep
endif
ifeq ($(MAKELEVEL), 2)
STEP = SecondStep
endif
ifeq ($(MAKELEVEL), 3)
STEP = ThirdStep
endif
ifeq ($(SRC),)
STEP = NullStep
endif
ifeq ($(INIT),init)
STEP = InitStep
endif

#-------------------  Rules for Handling compilation --------------------------#

# General rule. Switches execution based on each step.
all: $(STEP)

# Compile the program.

$(TARGET): $(OBJS)
	$(CC) $(CPPFLAGS) -o $@ $^

# Create symbolic links to source files.
$(SRC): % :
	rm -f $(LDIR)/$@-dir
	ln -s $(shell find $(SEDIR) -maxdepth 1 -name '$*.cpp' -print | \
	sed 's,\(.*/\)$*\.cpp,\1,g' | sed 's, ,$(sp),g' | \
	sed '\,\./, s,\./,../../,') '$(LDIR)/$@-dir'

# Compile object files from the source files. Also retrieves dependency info.
$(OBJS): $(ODIR)/%.o: %.cpp 
	$(CC) $(DEPFLAGS) '$(DDIR)/$*.d' '$(LDIR)/$*-dir/$(<F)'
	@echo '' >> '$(DDIR)/$*.d'
	@echo '$(LDIR)/$*-dir/$(<F):' >> '$(DDIR)/$*.d'
	$(CC) -c $(CPPFLAGS) '$(LDIR)/$*-dir/$(<F)' -o '$@'
	
# This line handles recompilation for dependency changes.
ifeq ($(STEP), SecondStep)
-include $(DEPS)
endif

#---------------------- Compilation Sequence ----------------------------------#

# Show source directories
ZerothStep: Tag0 
	($(MAKE))

# Make symbolic links
FirstStep: Tag1 $(SRC) 
	($(MAKE))

# Check for updates to source code and dependencies.
SecondStep: Tag2 $(TARGET)
	($(MAKE))

ThirdStep: Tag3

# No source files found.
NullStep: TagNull

# No source folder found.
InitStep: TagInit

#-------------------- Progress Information-------------------------------------#

Tag0:
	@echo  
	@echo -----------------------------------------------
	@echo . \*\*\*\* Preparing $(NAME) For Build \*\*\*\* .
	@echo -----------------------------------------------
	@echo  
	@echo Source directories: '$(SEDIR)'

Tag1:
	@echo  
	@echo -----------------------------------------------
	@echo . \*\*\*\* Getting Symbolic Links to Sources \*\*\*\* .
	@echo -----------------------------------------------
	@echo 

Tag2:
	@echo 
	@echo -----------------------------------------------
	@echo . \*\*\*\* Compiling and Linking Source Code \*\*\*\* .
	@echo -----------------------------------------------
	@echo 

Tag3:
	@echo 
	@echo -----------------------------------------------
	@echo . \*\*\*\* FINISHED \*\*\*\* .
	@echo -----------------------------------------------
	@echo
	
TagNull:
	@echo No source files found! \
	If you are not starting a new project, \
	did you place your source files in $(SDIR) \
	or put their path in the variable EXTDIRS?
	
TagInit:
	@echo No source folder found. Initializing folder hierarchy.

#-----------------Rules which do not make files--------------------------------#

.PHONY: clean cleaner $(STEPS) $(TAGS)

#-----------------Rules for cleaning up the directory--------------------------#

clean: 
	rm -f $(TARGET) $(OBJS) $(LNKPATH)
	
cleaner:
	rm -r $(ODIR)
	rm -f $(TARGET)
	
#------------------------------------------------------------------------------#