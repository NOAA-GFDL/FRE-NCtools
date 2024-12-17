#include <version.h>
#include <stdio.h>

/*
  Sinmple utility to print the FRE-NCtools version
  string.  The user must pass in the command name
*/

int main(int argc, char *argv[])
{
    if (argc > 1)
    {
        print_version(argv[1]);
    }
    return 0;
}