CFLAGS= -O3 -std=c++17

FILES= main.o graph.o io.o

hegnel_vc: $(FILES)
	g++ $(CFLAGS) -o hegnel_vc main.o graph.o io.o

main.o: main.cpp
	g++ $(CFLAGS) -c main.cpp -o main.o

graph.o: graph.cpp graph.h
	g++ $(CFLAGS) -c graph.cpp -o graph.o

io.o: io.cpp io.h
	g++ $(CFLAGS) -c io.cpp -o io.o

clean:
	rm -f hegnel_vc $(FILES)
