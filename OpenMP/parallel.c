#include <stdio.h>

int main() {

#pragma omp parallel
    {
	printf("Ahoj\n");
    }
    
    return 0;
}
