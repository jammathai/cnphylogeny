TESTS := $(basename $(wildcard *.c))

all: $(TESTS)

%: %.c ../cnphylogeny.c
	$(CC) -lm -g -o $@ $< ../cnphylogeny.c

.PHONY: clean
clean:
	rm -f $(TESTS)
