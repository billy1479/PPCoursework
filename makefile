CXX = g++
CXXFLAGS = -std=c++11 -Wall
SRC = main.cc
TARGET = runme
TESTVARIABLE1 = [1.0 0.0 0.0] [0.0 1.0 0.0] [0.0 0.0 1.0]

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)
	rm -f result.txt

test: 
	./$(TARGET) $(TESTVARIABLE1)
