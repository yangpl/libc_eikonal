.PHONY: clean all
CFLAS = -Wall -g -o2
TGT = main
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)
LIB = -lm

all: $(TGT)

$(TGT): $(OBJ)
	gcc $(CFLAGS) $(OBJ) -o $@ $(LIB)
	
%.o: %.c
	gcc $(CFLAGS) -c $< -o $@ $(LIB)
	
clean:
	rm -f *.o main *.bin *.eps *.txt *~
