TESTS := $(basename $(wildcard *.c))

all: $(TESTS)

%: %.c ../src/cnphylogeny.c
	$(CC) -lm -g -o $@ $< ../src/cnphylogeny.c

.PHONY: clean
clean:
	rm -f $(TESTS)
