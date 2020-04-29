all: basic_SPN.o linear_cryptanalysis_lib.o break-basic_SPN.o
	gcc -pthread basic_SPN.o linear_cryptanalysis_lib.o break-basic_SPN.o -o break-basic_SPN

spn: basic_SPN.c basic_SPN.h
	gcc -c basic_SPN.c

lib: linear_cryptanalysis_lib.c linear_cryptanalysis_lib.h
	gcc -c linear_cryptanalysis_lib.c

break: break-basic_SPN.c break-basic_SPN.h
	gcc -c break-basic_SPN.c

clean:
	rm *.o break-basic_SPN