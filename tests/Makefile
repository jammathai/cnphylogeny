TESTS := $(basename $(wildcard *.c))

all: $(TESTS)

%: %.c
	$(CC) -g -o $@ $<

.PHONY: clean
clean:
	rm $(TESTS)
