#include <stdlib.h>

double drand48_() {
   return drand48();
}

void srand48_(int *pseedval) {
   (void) srand48((long) *pseedval);
}