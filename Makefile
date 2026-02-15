
build:
	mkdir -p build
	g++ -O3 -march=native  src/config.cpp src/data.cpp src/estimation.cpp src/parameters.cpp src/util.cpp src/statistics.cpp src/cache.cpp src/dynamic.cpp src/simulator.cpp src/main.cpp -o build/scache

install: build
	cp build/scache ./scache
	mkdir -p output

clean:
	rm -r ./build ./output scache
