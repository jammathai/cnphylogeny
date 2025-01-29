OBJS := $(patsubst src/%.c, build/%.o, $(wildcard src/*.c))

build/cnphylogeny: $(OBJS)
	$(CC) -lm -o $@ $^

build/%.o: src/%.c | build
	$(CC) -c -o $@ $<

build:
	mkdir build

.PHONY: clean
clean:
	rm -rf build
