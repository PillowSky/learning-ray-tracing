all:
	g++ -std=c++14 -fopenmp -Wall -march=native -O3 main.cpp `pkg-config --cflags --libs Magick++` -o main.out

clean:
	rm main.out