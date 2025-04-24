# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -I /users/cnwosu2/ParallelComp/rapidjson

# Linker flags
LDFLAGS = -lcurl

# Output file
OUTPUT = graph_crawler

# Source files
SRC = GraphCrawler.cpp

# Default target
all: $(OUTPUT)

# Compile and link
$(OUTPUT): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(SRC) $(LDFLAGS)

# Clean up compiled files
clean:
	rm -f $(OUTPUT)
