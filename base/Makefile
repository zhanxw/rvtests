
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
   HOME = /home/zhanxw
endif
ifeq ($(UNAME), Darwin)
   HOME = /Users/zhanxw
endif

EXEC = lib-base.a
BASE = Argument Exception IO OrderedMap Regex TypeConversion Utils
OBJ = $(BASE:=.o)

DEFAULT_CXXFLAGS = $(STATGEN_INC)

TARGET: CXXFLAGS = $(DEFAULT_CXXFLAGS)
TARGET: $(EXEC)

all: $(TARGET)

# 'make release' will: clean current build; build with -O4 flags; backup executable files to release/ folder
# NOTE: don't use -j flag!
release: CXXFLAGS = -O4 $(DEFAULT_CXXFLAGS)
release: clean $(EXEC) backup_release

debug: CXXFLAGS = -g $(DEFAULT_CXXFLAGS)
debug: $(EXEC)

$(EXEC): $(OBJ)
	ar rcs $(EXEC) $(OBJ)

%: %.cpp # remove implicit Makefile rule of converting ABC.cpp to ABC
%.o: %.cpp %.h
	g++ -c $(CXXFLAGS) $<

clean:
	rm -rf $(EXEC) *.o
