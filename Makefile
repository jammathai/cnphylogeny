OBJS := $(patsubst src/%.c, build/%.o, $(wildcard src/*.c))

build/cnphylogeny: $(OBJS)
	$(CC) -O5 -lm -o $@ $^

build/%.o: src/%.c | build
	$(CC) -O5 -c -o $@ $<

build:
	mkdir build

.PHONY: clean
clean:
	rm -rf build
