all:
	g++ hw1.cpp -o hw1 -std=c++17
	g++ hw1-o.cpp -o hw1-o -std=c++17
debug:
	g++ hw1.cpp -o hw1 -std=c++17 -DDEBUG
	g++ hw1-o.cpp -o hw1-o -std=c++17 -DDEBUG
clean:
	rm hw1 hw1-o