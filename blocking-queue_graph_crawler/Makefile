CXX= g++
CXXFLAGS= -std=c++17 -Wall -I/usr/local/include
LDFLAGS= -lcurl -pthread

TARGET= main
SRCS= main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f main *.o output/*.txt
