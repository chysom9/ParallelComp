CXXFLAGS = -std=c++17 -I /Users/chi/Desktop/Spring2025/ParallelComp/rapidjson/include
LDFLAGS=-lcurl
LD=g++
CC=g++

all: level_client par_level_client

level_client: level_client.cpp
	$(LD) $(CXXFLAGS) level_client.cpp -o level_client $(LDFLAGS)

par_level_client: par_level_client.cpp
	$(LD) $(CXXFLAGS) par_level_client.cpp -o par_level_client $(LDFLAGS)



clean:
	-rm -f level_client par_level_client
