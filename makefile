CC=gcc

INCLUDE_FLAGS= -I. -I/usr/include/python3.10/ -lm -lpython3.10
WARNING_FLAGS=-Wall -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations
CFLAGS = $(INCLUDE_FLAGS) $(WARNING_FLAGS)

all: test linux 

test:
	$(CC) ./smile.c -o smile.out $(CFLAGS)

linux:
	$(CC) -c -fPIC ./smile.c -o smile.o $(CFLAGS)
	$(CC) smile.o -shared -o smile.so

clean: 
	rm ./*.out
	rm ./*.so
	rm ./*.o