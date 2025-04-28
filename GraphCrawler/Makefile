# Compiler and flags
CXX = g++
<<<<<<< HEAD
CXXFLAGS = -std=c++17 -I /u sers/cnwosu2/ParallelComp/rapidjson/include
=======
CXXFLAGS = -std=c++17 -I /users/cnwosu2/ParallelComp/rapidjson
>>>>>>> 58a9893717612a279998a43bc05edb23821e4b9a

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
