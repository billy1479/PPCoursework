CXX = g++
CXXFLAGS = -std=c++11 -Wall
SRC = main.cc
TARGET = runme
TESTVARIABLE1 = "[1,0,0]"
TESTVARIABLE2 = "[0,1,0]"
TESTVARIABLE3 = "[0,0,1]"

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)

test: 
	./$(TARGET) $(TESTVARIABLE1) $(TESTVARIABLE2) $(TESTVARIABLE3)
