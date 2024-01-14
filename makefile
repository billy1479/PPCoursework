CXX = g++
CXXFLAGS = -std=c++11 -Wall
SRC = main.cc
TARGET = runme
TEST = [15.0 -7.0 -7.0] [-7.0 15.0 -7.0] [-7.0 -7.0 15.0]

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)
	rm -f result.txt

test: 
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)
	./$(TARGET) $(TEST)
