
# 1. Put this file in the same folder as your 'driver' code
#    (the code containing the 'main' function).

# 2. Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it)
#    DMT_DIR is the directory with the DMT header files.
LIBRARY_DIR=/n/home06/jkemp/ITensor
DMT_DIR=/n/home06/jkemp/SharedDMT

# 3. If your 'main' function is in a file called 'myappname.cc', then
#    set APP to 'myappname'. Running 'make' will compile the app.
#    Running 'make debug' will make a program called 'myappname-g'
#    which includes debugging symbols and can be used in gdb (Gnu debugger);
APP=xyz

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
HEADERS=

# 5. For any additional .cc files making up your project,
#    add their full filenames here.
CCFILES=$(APP).cc

#################################################################
#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

CCFLAGS += -Wno-unused-result -Wno-deprecated-declarations -I$(DMT_DIR)
CCGFLAGS += -I$(DMT_DIR)



TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h

#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

ROBJECTS=$(patsubst %.cc,%.o, $(RCCFILES))
RGOBJECTS=$(patsubst %,.debug_objs/%, $(ROBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -Wno-ignored-attributes -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) -Wno-ignored-attributes $(OBJECTS) -o $(APP) $(LIBFLAGS)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)


clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g

run:
	./$(APP)

rund:
	./$(APP)-g

brun:
	$(APP)
	./$(APP)

mkdebugdir:
	mkdir -p .debug_objs
