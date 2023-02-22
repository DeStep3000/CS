#include <stdio.h>

int main() {
    printf("Hello, what's your name?\n");
    char name[50];
    scanf("%s", name);
    printf("Hello, %s. Welcome to CS (Cubic_Spline).\n"
           "This program is designed to find cubic splines and their intersections with other splines.\n"
           "Let's enter the number of our points:\n", name);
    int n;
    scanf("%d", &n);
    while (n <= 1){
        printf("Please, enter a positive integer starting from 2:\n");
        scanf("%d", &n);
    }
    return 0;
}
