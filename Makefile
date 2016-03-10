TARGET = example
OBJS = example.o
ifeq ($(OS), Windows_NT)
	EXEEXT = .exe
endif
CPPFLAGS = -std=c++11
LDFLAGS = -lstdc++ -lm

$(TARGET)$(EXEEXT): $(OBJS)
