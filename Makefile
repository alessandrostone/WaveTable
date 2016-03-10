TARGET = example.exe
OBJS = example.o
CPPFLAGS = -std=c++11
LDFLAGS = -lstdc++ -lm

$(TARGET): $(OBJS)
