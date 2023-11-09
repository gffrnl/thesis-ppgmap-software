header_files = $(wildcard ./*.h) \
               $(wildcard ../include/*.h) \
               $(wildcard ../include/frlap/*.h) \
               $(wildcard ../src/*.h) \
               $(wildcard ../src/frlap1dq/*.h) \
               ../src/toeplitz/toeplitz.h

src_files = $(wildcard ../src/*.c) \
            $(wildcard ../src/frlap1dq/*.c) \
            ../src/toeplitz/toeplitz.c

obj_files = $(src_files:.c=.o)

CFLAGS += -g
CPPFLAGS += -I../include -I../src/
LDLIBS += -lc -lm -lfftw3

all: bench

run: bench
	./bench

#bench__naive__: $(header_files) $(obj_files) bench__naive__.o
#	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $(obj_files) bench__naive__.o

bench: $(header_files) $(obj_files) bench.o
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $(obj_files) bench.o

#bencht: $(header_files) $(obj_files) bencht.o
#	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $(obj_files) bencht.o


%.o: %.c $(header_files)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

#frlap/toeplitz/toeplitz.o: frlap/toeplitz/toeplitz.c frlap/toeplitz/toeplitz.h

clean:
#	rm -f bencht
	rm -f bench
	rm -f *.o
	rm -f ../src/*.o
	rm -f ../src/frlap1dq/*.o
	rm -f ../src/toeplitz/toeplitz.o

purge: clean
	rm -f *~
#	rm -f ../frlap/*~
#	rm -f frlap/frlap1dq/*~


print_headers:
	ls -la $(header_files)

print_sources:
	ls -la $(src_files)

.PHONY: all clean purge print_sources print_headers
