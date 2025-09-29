all:
	g++ hw1.cpp -o hw1 -std=c++17
debug:
	g++ hw1.cpp -o hw1 -std=c++17 -DDEBUG
clean:
	rm hw1 hw1-o