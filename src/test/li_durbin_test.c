#include "unittests.h"
#include <testlib.h>

int main(void)
{
    TL_BEGIN_TEST_SUITE("Li & Durbin test");
#ifdef GEN_UNIT_TESTS // unit testing of static functions...
    // Unit tests
    TL_RUN_TEST(build_ld_tables);
    TL_RUN_TEST(ld_iterator);
#endif

    TL_END_SUITE();
}
