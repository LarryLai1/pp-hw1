all:
	g++ hw1.cpp -o hw1 -std=c++17
debug:
	g++ hw1.cpp -o hw1 -std=c++17 -DDEBUG -g
clean:
	rm hw1 hw1-o