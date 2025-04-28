#CXXFLAGS=-I ~/prgs/rapidjson/include
LDFLAGS=-lcurl
LD=g++
CC=g++

all: client

client: client.o
	$(LD) $< -o $@ $(LDFLAGS)

clean:
	-rm client client.o
