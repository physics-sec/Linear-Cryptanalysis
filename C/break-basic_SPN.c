#include <stdio.h>
#include "basic_SPN.h"
#include "linear_cryptanalysis_lib.h"

int main()
{
    struct state** linear_aproximations = analize_cipher();
    for (int i = 0; i < 30; i++)
    {
        printState(*linear_aproximations[i]);
    }
    struct state state = *linear_aproximations[0];
    int result = get_xor(0, 29594, 0, state);
    printf("result: %d", result);
    freeMem(linear_aproximations);
    return 0;
}
