#include <stdio.h>
#include <sys/time.h>

// Factorial function of int n without recursion
unsigned long long int factorial(int n)
{
    unsigned long long int result = 1;
    for (int i = 1; i <= n; i++)
    {
        result *= i;
    }
    return result;
}

int main()
{
    int n = 20;

    // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);

    // For loop a 10000 times to get a more accurate time
    for (int i = 0; i < 10000; i++)
    {
        factorial(n);
    }

    // Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = (seconds + microseconds*1e-6)/10000;

    printf("Factorial of %d = %llu\n", n, factorial(n));
    printf("Time measured: %.9f seconds.\n", elapsed);

    return 0;
}