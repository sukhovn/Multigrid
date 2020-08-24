.PHONY: all clean

PROGRAM = multigrid
MAIN = test
CORE_DIR = Core

CORE_SRC := $(wildcard $(CORE_DIR)/*.cpp)
CORE_BINS := $(patsubst $(CORE_DIR)/%.cpp, obj/core/%.o, $(CORE_SRC))

HEAD_FOLDERS = $(CORE_DIR)
HEAD_LINK = $(patsubst %, -I%, $(HEAD_FOLDERS))

CC = g++ # compiler to use
CFLAGS = -c -O3 $(HEAD_LINK) # flags to use at the compliation
LINKERFLAG = -lm

all: prog

obj:
	@echo "Making directory for binaries..."
	mkdir obj
	mkdir obj/core

obj/core/%.o: $(CORE_DIR)/%.cpp obj
	@echo "Compiling core binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/$(MAIN).o: $(MAIN).cpp obj
	@echo "Compiling the runner..."
	${CC} ${CFLAGS} -o obj/$(MAIN).o $(MAIN).cpp

prog: $(CORE_BINS) obj/$(MAIN).o
	@echo "Linking the tester..."
	${CC} -o $(PROGRAM) $(CORE_BINS) obj/$(MAIN).o $(LINKERFLAG)

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -vf $(PROGRAM)