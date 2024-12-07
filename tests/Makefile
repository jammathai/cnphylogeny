SRCS := $(wildcard *.c)
TESTS := $(basename $(SRCS))
FLAGS := -lm -g

all: $(TESTS)

%: %.c ../cnphylogeny.c
	$(CC) $(FLAGS) -o $@ $< ../cnphylogeny.c

clean:
	rm $(TESTS)
