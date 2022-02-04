#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_usage(const char *progname)
{
    fprintf(stderr,
            "Usage: %s -p genome\n       %s -d dist genome reads\n",
            progname, progname);
    exit(1);
}

int main(int argc, char const *argv[])
{
    if ((argc == 3) && strcmp("-p", argv[1]) == 0)
    {
        // preprocessing
        printf("Preprocessing genome %s.\n", argv[2]);
    }
    else if ((argc == 5) && strcmp("-d", argv[1]) == 0)
    {
        int d = atoi(argv[2]); // you can do better than atoi, but I'm feeling lazy
        const char *genome = argv[3];
        const char *reads = argv[4];
        printf("Mapping in genome %s for reads in %s within distance %d.\n",
               genome, reads, d);
    }
    else
    {
        print_usage(argv[0]);
    }

    return 0;
}
