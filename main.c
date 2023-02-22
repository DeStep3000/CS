#include <stdio.h>

void hello(){
    //  Делаем приветсвие
    printf("Hello, what's your name?\n");
    char name[50];
    scanf("%s", name);
    printf("Hello, %s. Welcome to CS (Cubic_Spline).\n"
           "This program is designed to find cubic splines and their intersections with other splines.\n"
           "Let's enter the number of our points:\n", name);

//  Вводим количество точек
    double n0;
    scanf("%lf", &n0);

//  Провереям, чтобы количество точек было целым числом и чтобы было как минимум 2 точки
    while (n0 <= 1 || n0 != (int) n0){
        printf("Please, enter a positive integer starting from 2:\n");
        scanf("%lf", &n0);
    }
    int n = (int) n0;

//   Вводим x
    printf("Enter the function arguments(x):\n");
    double x[n];
    for (int i = 0; i < n; i++){
        scanf("%lf", &x[i]);
    }

//  Вводим y
    printf("Enter the function values(y):\n");
    double y[n];
    for (int i = 0; i < n; i++){
        scanf("%lf", &y[i]);
    }
}

int main() {
    hello();
    return 0;
}
