#include <stdio.h>
#include "ran2.c"

long seed = -1;

int main() {
	int i;

	for(i=0; i<10; i++) {
		printf("%f\n", ran2(&seed));
	}
}
