CXX = g++
CXXFLAGS = -std=c++17 -Wall -pthread -I/rapidjson/include

all: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o exec -lcurl

clean:
	rm -f exec *.o output/*.txt