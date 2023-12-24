CXX = g++
CXXFLAGS = -std=c++11 -Wall
SRC = main.cc
TARGET = runme

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)