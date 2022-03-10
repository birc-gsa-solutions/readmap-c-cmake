#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "unittests.h"
#include <testlib.h>

int main(void)
{
    TL_BEGIN_TEST_SUITE("skew test");
#ifdef GEN_UNIT_TESTS // unit testing of static functions...
    TL_RUN_TEST(skew_test_len_calc);
#endif
    TL_END_SUITE();
}
