#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Константа, определяющая точность вычислений
#define EPSILON 1e-6
#define M_PI 3.141592653589793


// Структура, описывающая кубический сплайн
typedef struct {
    double a, b, c, d; // Коэффициенты для каждого интервала
    double x; // Левый конец интервала
} CubicSpline;


struct {
    double first;
    double second;
} DeStep[3] = {EPSILON, EPSILON, EPSILON, EPSILON, EPSILON, EPSILON};


// Функция, которая вычисляет коэффициенты для кубического сплайна
// на основе заданных координат x и y точек
void compute_spline_coefficients(double *x, double *y, int n, CubicSpline *spline) {
    double *h = malloc(n * sizeof(double));
    double *alpha = malloc(n * sizeof(double));
    double *l = malloc(n * sizeof(double));
    double *mu = malloc(n * sizeof(double));
    double *z = malloc(n * sizeof(double));
    int i;

    for (i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }

    for (i = 1; i < n - 1; i++) {
        alpha[i] = 3.0 / h[i] * (y[i + 1] - y[i]) - 3.0 / h[i - 1] * (y[i] - y[i - 1]);
    }

    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1;
    z[n - 1] = 0;
    spline[n - 1].c = 0;


    for (i = n - 2; i >= 0; i--) {
        spline[i].c = z[i] - mu[i] * spline[i + 1].c;
        spline[i].b = (y[i + 1] - y[i]) / h[i] - h[i] * (spline[i + 1].c + 2.0 * spline[i].c) / 3.0;
        spline[i].d = (spline[i + 1].c - spline[i].c) / (3.0 * h[i]);
        spline[i].a = y[i];
        spline[i].x = x[i];
    }

    free(h);
    free(alpha);
    free(l);
    free(mu);
    free(z);
}

// Функция, которая вычисляет расстояние между двумя точками
double compute_distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double aNDREW(double x) {
    if (x < 0)
        return -pow(-x, 1.0 / 3.0);
    return pow(x, 1.0 / 3.0);
}

void square(double a, double b, double c) {
    double discr = b * b - 4 * a * c;

    if (discr < 0) {
        for (int i = 0; i < 3; i++)
            DeStep[i].second = EPSILON;
    } else if (discr == 0) {
        DeStep[0].first = ((-b + sqrt(discr)) / 2 * a);
        DeStep[0].second = 0;
        for (int i = 1; i < 3; i++) { DeStep[i].second = EPSILON; }
    } else {
        printf("x = %f and x = %f \n", ((-b + sqrt(discr)) / 2 * a), ((-b - sqrt(discr)) / 2 * a));
        for (int i = 0; i < 2; i++) {
            DeStep[i].second = 0;
            DeStep[i].first = ((-b + pow((-1), i) * sqrt(discr)) / 2 * a);
        }
        DeStep[2].second = EPSILON;
    }
}

double max(double a, double b) {
    return (a > b ? a : b);
}

double min(double a, double b) {
    return (a < b ? a : b);
}

// Функция, которая проверяет, пересекаются ли два кубических сплайна, и если да, находит точки пересечения
int find_intersection(CubicSpline *spline1, CubicSpline *spline2, int n, int m, double *x1, double *x2,
                       double *y1, double *y2) {
    double right, left;
    for (int i = 0; i < n; i++) {
        for (int j = 0; i < m; i++) {
            right = max(x1[i], x2[j]);
            left = min(x1[i + 1], x2[j + 1]);
            if (left < right) {
                double A = spline1[i].d - spline2[j].d;
                double B = spline1[i].c - 3 * spline1[i].d * x1[i] - spline2[j].c + 3 * spline2[j].d * x2[j];
                double C = spline1[i].b - 2 * spline1[i].c * x1[i] + 3 * spline1[i].d * x1[i] * x1[i] - spline2[j].b +
                           2 * spline2[j].c * x2[j] - 3 * spline2[j].d * x2[j] * x2[j];
                double D = y1[i] - spline1[i].b * x1[i] + spline1[i].c * x1[i] * x1[i] -
                           spline1[i].d * x1[i] * x1[i] * x1[i] - y2[j] + spline2[j].b * x2[j] -
                           spline2[j].c * x2[j] * x2[j] + spline2[j].d * x2[j] * x2[j] * x2[j];

                if (A == 0 && B != 0)
                    square(B, C, D);

                else if (A == 0 && B == 0) {// Cx + D = 0
                    DeStep[0].first = (-D / C);
                    DeStep[0].second = 0;
                } else {
                    double p = (3.0 * A * C - B * B) / (3.0 * A * A);
                    double q = (2.0 * B * B * B - 9.0 * A * B * C + 27.0 * A * A * D) / (27.0 * A * A * A);
                    double S = (q * q / 4.0) + (p * p * p / 27.0);

                    double F;
                    if (q == 0)
                        F = M_PI / 2.0;
                    if (q < 0)
                        F = atan(-2.0 * sqrt(-S) / q);
                    if (q > 0)
                        F = atan(-2.0 * sqrt(-S) / q) + M_PI;

                    for (int i = 0; i < 3; i++)
                        DeStep[i].first = DeStep[i].second = 0;

                    if (S < 0) {
                        DeStep[0].first = 2.0 * sqrt(-p / 3.0) * cos(F / 3.0) - B / (3.0 * A);
                        DeStep[1].first = 2.0 * sqrt(-p / 3.0) * cos((F / 3.0) + 2.0 * M_PI / 3.0) - B / (3.0 * A);
                        DeStep[2].first = 2.0 * sqrt(-p / 3.0) * cos((F / 3.0) + 4.0 * M_PI / 3.0) - B / (3.0 * A);
                    }

                    if (S == 0) {
                        DeStep[0].first = 2.0 * aNDREW(-q / 2.0) - B / (3.0 * A);
                        DeStep[1].first = -aNDREW(-q / 2.0) - B / (3.0 * A);
                        DeStep[2].first = -aNDREW(-q / 2.0) - B / (3.0 * A);
                    }

                    if (S > 0) {
                        double temp1 = aNDREW((-q / 2.0) + sqrt(S)) + aNDREW((-q / 2.0) - sqrt(S));
                        double temp2 = aNDREW((-q / 2.0) + sqrt(S)) - aNDREW((-q / 2.0) - sqrt(S));
                        DeStep[0].first = temp1 - B / (3.0 * A);
                        DeStep[1].first = -temp1 / 2.0 - B / (3.0 * A);
                        DeStep[1].second = sqrt(3) * temp2 / 2.0;
                        DeStep[2].first = -temp1 / 2.0 - B / (3.0 * A);
                        DeStep[2].second = -sqrt(3) * temp2 / 2.0;
                    }
                }
            }
        }
    }

    int flag = 0;
    int num = 1;

    // Если сплайны пересекаются, выводим координаты точек пересечения
    printf("Spline_1 and spline_2 intersect at:\n");
    for (int i = 0; i < 3; i++) {
        if (DeStep[i].second == 0 && DeStep[i].first != EPSILON) {
            printf("%d.  x = %f\n", num, DeStep[i].first);
            flag = 1;
            num++;
        }
    }
    return flag;
}

// Функция, которая находит минимальное расстояние между двумя кубическими сплайнами
double find_min_distance(CubicSpline *spline1, CubicSpline *spline2, int n, int m) {
    double min_distance = INFINITY;
    int i, j;
    // Вычисляем расстояние между каждой парой точек на обоих сплайнах
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < m - 1; j++) {
            double x1 = spline1[i].x;
            double y1 = spline1[i].a + spline1[i].b * fabs((x1 - spline1[i + 1].x)) +
                        spline1[i].c * pow(x1 - spline1[i + 1].x, 2) +
                        spline1[i].d * fabs(pow(x1 - spline1[i].x, 3));
            double x2 = spline2[j].x;
            double y2 = spline2[j].a + spline2[j].b * fabs((x2 - spline2[j + 1].x)) +
                        spline2[j].c * pow(x2 - spline2[j + 1].x, 2) +
                        spline2[j].d * fabs(pow(x2 - spline2[j + 1].x, 3));
            double distance = compute_distance(x1, y1, x2, y2);
            if (distance < min_distance) {
                min_distance = distance;
            }
        }
    }

    return min_distance;
}

// Совместимость сплайнов
int is_cubic_spline_same(CubicSpline *spline1, CubicSpline *spline2, int n, int m) {
    if (n != m) {
        return 0;
    }
    for (int i = 0; i < n - 1; i++) {
        if (fabs(spline1[i].a - spline2[i].a) > EPSILON ||
            fabs(spline1[i].b - spline2[i].b) > EPSILON ||
            fabs(spline1[i].c - spline2[i].c) > EPSILON ||
            fabs(spline1[i].d - spline2[i].d) > EPSILON) {
            return 0;
        }
    }
    return 1;
}

// Функция для вычисления значения кубического сплайна в точке x
double evaluate_cubic_spline(double x, double *splX, CubicSpline *spline, int n) {
    if (x <= splX[0])
        return spline[0].a + spline[0].b * (x - splX[0]) + spline[0].c * pow(x - splX[0], 2) +
               spline[0].d * pow(x - splX[0], 3);

    if (x >= splX[n - 2])
        return spline[n - 2].a + spline[n - 2].b * (x - splX[n - 2]) + spline[n - 2].c * pow(x - splX[n - 2], 2) +
               spline[n - 2].d * pow(x - splX[n - 2], 3);

    int i;
    for (i = 0; i < n - 1; i++) {
        if (x >= splX[i] && x < splX[i + 1])
            break;
    }


    double h = x - splX[i];
    double y = spline[i].a + spline[i].b * h + spline[i].c * pow(h, 2) + spline[i].d * pow(h, 3);
    return y;
}


void sort(double *arr_x, double *arr_y, int n, int k) {
    int noSwap;
    int tmp;
    for (int i = n - 1; i >= 0; i--) {
        noSwap = 1;
        for (int j = 0; j < i; j++) {
            if (arr_x[j] > arr_x[j + 1]) {
                tmp = arr_x[j];
                arr_x[j] = arr_x[j + 1];
                arr_x[j + 1] = tmp;
                noSwap = 0;
                tmp = arr_y[j];
                arr_y[j] = arr_y[j + 1];
                arr_y[j + 1] = tmp;
            }
        }
        if (noSwap == 1)
            break;
    }
    printf("Spline_%d sorted value of x and y:\n", k);
    for (int i = 0; i < n; i++) {
        printf("         %lf %lf\n", arr_x[i], arr_y[i]);
    }
    printf("------------------------------------------------------------------------------------------\n");

}

int main() {
    // ЕСЛИ ЭТО ВСЁ РАСКОМЕНТИРОВАТЬ, КОД ВЫВОДИТ НЕ "Spline 1 and spline 2 are the same."
    // ЭТО ВСЁ ИДЕТ, КОГДА МЫ МЕНЯЕМ ПО ДРУГОМУ ВВОДИМ ЧИСЛА В МАССИВ

    //  Делаем приветсвие
    printf("Hello, what's your name?\n");
    char name[50];
    scanf("%s", name);
    printf("Hello, %s. Welcome to CS (Cubic_Spline).\n"
           "This program is designed to find cubic splines and their intersections with other splines.\n"
           "------------------------------------------------------------------------------------------\n"
           "Let's enter the number of our points to splain_1:\n", name);

    //  Вводим количество точек для 1 сплайна
    double n0;
    scanf("%lf", &n0);

    //  Провереям, чтобы количество точек было целым числом и чтобы было как минимум 2 точки
    while (n0 <= 1 || n0 != (int) n0) {
        printf("Please, enter a positive integer starting from 2:\n");
        scanf("%lf", &n0);
    }
    int n = (int) n0;

    //   Вводим x
    printf("Enter the function arguments(x):\n");
    double x1[n];
    for (int i = 0; i < n; i++) {
        scanf("%lf", &x1[i]);
    }

    //  Вводим y
    printf("Enter the function values(y):\n");
    double y1[n];
    for (int i = 0; i < n; i++) {
        scanf("%lf", &y1[i]);
    }
    printf("------------------------------------------------------------------------------------------\n");

    // Также делаем и для второго сплайна
    printf("Let's enter the number of our points to splain_2:\n");

    double m0;
    scanf("%lf", &m0);

    while (m0 <= 1 || m0 != (int) m0) {
        printf("Please, enter a positive integer starting from 2:\n");
        scanf("%lf", &m0);
    }
    int m = (int) m0;

    printf("Enter the function arguments(x):\n");
    double x2[m];
    for (int i = 0; i < m; i++) {
        scanf("%lf", &x2[i]);
    }

    printf("Enter the function values(y):\n");
    double y2[m];
    for (int i = 0; i < m; i++) {
        scanf("%lf", &y2[i]);
    }

//    int n = 5;
//     // инициализация координат x и y точек сплайнов
//    double x1[] = {0, 1, 2, 3, 4};
//    double y1[] = {0, 1, 4, 9, 16};
//
//    int m = 5;
//    double x2[] = {0, 1, 2, 3, 4};
//    double y2[] = {0, 1, 4, 9, 16};
    printf("------------------------------------------------------------------------------------------\n");
    // сортировка сплайнов
    sort(x1, y1, n, 1);
    sort(x2, y2, m, 2);

    // создание сплайнов
    CubicSpline spline1[n];
    CubicSpline spline2[m];

    // вычисление коэффициентов
    compute_spline_coefficients(x1, y1, n, spline1);
    compute_spline_coefficients(x2, y2, m, spline2);

    //Выводим коэффициенты
    for (int i = 0; i < n - 1; i++) {
        printf("Spline_1[%d]: a = %lf, b = %lf, c = %lf, d = %lf, x = %lf\n", i, spline1[i].a,
               spline1[i].b, spline1[i].c, spline1[i].d, spline1[i].x);
    }
    printf("\n");
    for (int i = 0; i < m - 1; i++) {
        printf("Spline_2[%d]: a = %lf, b = %lf, c = %lf, d = %lf, x = %lf\n", i, spline2[i].a,
               spline2[i].b, spline2[i].c, spline2[i].d, spline2[i].x);
    }
    printf("------------------------------------------------------------------------------------------\n");

    // Спрашиваем у пользователя, в каком сплайне мы хотим вычислить значение точки сплайна
    // И вычисляем
    double new_x;
    int spline_num;
    double new_y;
    printf("Enter the spline number, which you would like to work on:\n");
    scanf("%d", &spline_num);

    while (spline_num != 1 && spline_num != 2) {
        printf("The spline number can be only 1 or 2:\n");
        scanf("%d", &spline_num);
    }

    printf("Enter the point, which function value you want to know:\n");
    scanf("%lf", &new_x);

    if (spline_num == 1) {
        new_y = evaluate_cubic_spline(new_x, x1, spline1, n);
    } else {
        new_y = evaluate_cubic_spline(new_x, x2, spline2, m);
    }
    printf("\n");
    printf("The value of function in point %lf = %lf\n", new_x, new_y);
    printf("------------------------------------------------------------------------------------------\n");

    // проверка совпадения сплайнов
    int is_same_spline = is_cubic_spline_same(spline1, spline2, n, m);
    if (is_same_spline) {
        printf("Spline_1 and spline_2 are the same.\n");
        return 0;
    }

    // проверка на пересечение сплайнов
    // ВОТ НАД ЭТИМ НАДО РАБОТАТЬ
    int intersect = find_intersection(spline1, spline2, n, m, x1, x2, y1, y2);

    if (intersect == 0) {
        // Если сплайны не пересекаются, находим минимальное расстояние между ними
        double min_distance = find_min_distance(spline1, spline2, n, m);
        printf("The minimum distance between spline 1 and spline 2 is %f.\n", min_distance);
    }
    return 0;
}
