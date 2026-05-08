#include <stdio.h>

// Something to bind in libreltrans so that the linker puts it all together.
extern double disco_(double *);

int main() {
    double spin = 0.998;
    double isco = disco_(&spin);
    printf("Success: %lf\n", isco);
    return 0;
}
