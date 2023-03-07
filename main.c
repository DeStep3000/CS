#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Константа, определяющая точность вычислений
#define EPSILON 1e-6
#define PI 3.141592653589793
#define PI_2 (2 * PI)

// Структура, описывающая кубический сплайн
typedef struct {
    double a, b, c, d; // Коэффициенты для каждого интервала
    double x; // Левый конец интервала
} CubicSpline;

//Структура, описывающая координаты точек
typedef struct {
    double x;
    double y;
} Coords;

//Собственные функции min и max соответственно
double min(double a, double b) {
    return (a < b ? a : b);
}

double max(double a, double b) {
    return (a > b ? a : b);
}

// Функция, которая вычисляет коэффициенты для кубического сплайна
// на основе заданных координат x и y точек
void compute_spline_coefficients(double *x, double *y, int n, CubicSpline *spline) {
    // Выделяем память
    double *h = malloc(n * sizeof(double));
    double *alpha = malloc(n * sizeof(double));
    double *l = malloc(n * sizeof(double));
    double *mu = malloc(n * sizeof(double));
    double *z = malloc(n * sizeof(double));
    int i;

    // Находим приращение аргумента
    for (i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }

    // Находим alpha - вектор правой части системы линейных уравнений
    for (i = 1; i < n - 1; i++) {
        alpha[i] = 3.0 / h[i] * (y[i + 1] - y[i]) - 3.0 / h[i - 1] * (y[i] - y[i - 1]);
    }
    // Начальные значения
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;
    // Алгоритм прогонки
    for (i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1;
    z[n - 1] = 0;
    spline[n - 1].c = 0;

    // Заполняем значения коэффициентов сплайна
    for (i = n - 2; i >= 0; i--) {
        spline[i].c = z[i] - mu[i] * spline[i + 1].c; // С содержит решение СЛАУ
        spline[i].b = (y[i + 1] - y[i]) / h[i] - h[i] * (spline[i + 1].c + 2.0 * spline[i].c) / 3.0; // Вычисляется по формуле из первой производной
        spline[i].d = (spline[i + 1].c - spline[i].c) / (3.0 * h[i]); // Вычисляется по формуле из третьей производной
        spline[i].a = y[i]; // равно значению функции в точке i
        spline[i].x = x[i];
    }

    // Освобождаем память
    free(h);
    free(alpha);
    free(l);
    free(mu);
    free(z);
}

// Функция для вычисления значения кубического сплайна в точке x
// Нужна для вывода
double evaluate_cubic_spline(double x, double *splX, CubicSpline *spline, int n) {
    // Вычисление точки, находящейся левее нашего сплайна
    if (x < splX[0])
        return spline[0].a + spline[0].b * (x - splX[0]) + spline[0].c * pow(x - splX[0], 2) +
               spline[0].d * pow(x - splX[0], 3);

    // Вычисление точки, находящейся правее нашего сплайна
    if (x > splX[n - 1])
        return spline[n - 2].a + spline[n - 2].b * (x - splX[n - 1]) + spline[n - 2].c * pow(x - splX[n - 1], 2) +
               spline[n - 2].d * pow(x - splX[n - 1], 3);

    // Находим в каком промежутке находится наша точка x
    int i;
    for (i = 0; i < n - 1; i++) {
        if (x >= splX[i] && x <= splX[i + 1])
            break;
    }

    // Вычисляем значение сплайна в точке x
    double h = x - splX[i];
    double y = spline[i].a + spline[i].b * h + spline[i].c * pow(h, 2) + spline[i].d * pow(h, 3);
    return y;
}

// Функция для вычисления значения кубического сплайна в точке x на некотором сплайне с известной x0
// Нужна для вычислений в программе
double calculate_point(CubicSpline *coefs, double x, double x0) {
    return coefs->a + coefs->b * (x - x0) + coefs->c * pow(x - x0, 2) + coefs->d * pow(x - x0, 3);
}

// Функция для решения квадратного уравнения
int solve_quadratic(double *x, double a, double b, double c) {
    double d;
    d = b * b - 4. * a * c;
    if (d == 0.) {
        x[0] = x[1] = -b / (2. * a);
        return 1;
    }
    if (d < 0.) {
        return 0;
    }
    // 2 корня: избегаем вычитания 2 близких чисел
    if (b >= 0)
        d = (-0.5) * (b + sqrt(d));
    else
        d = (-0.5) * (b - sqrt(d));
    x[0] = d / a;
    x[1] = c / d;
    return 2;
}

// Функция для решения кубического уравнения методом Кардано
int solve_cubic(double *x, double a, double b, double c) {
    double q, r, r2, q3;
    q = (a * a - 3. * b) / 9.;
    r = (a * (2. * a * a - 9. * b) + 27. * c) / 54.;
    r2 = r * r;
    q3 = q * q * q;
    if (r2 < q3) {
        double t = acos(r / sqrt(q3));
        a /= 3.;
        q = -2. * sqrt(q);
        x[0] = q * cos(t / 3.) - a;
        x[1] = q * cos((t + PI_2) / 3.) - a;
        x[2] = q * cos((t - PI_2) / 3.) - a;
        return (3);
    } else {
        double aa, bb;
        if (r <= 0)
            r = -r;
        aa = -pow(r + sqrt(r2 - q3), 1. / 3.);
        if (aa != 0)
            bb = q / aa;
        else
            bb = 0.;
        a /= 3.;
        q = aa + bb;
        r = aa - bb;
        x[0] = q - a;
        x[1] = (-0.5) * q - a;
        x[2] = (sqrt(3.) * 0.5) * fabs(r);
        if (x[2] == 0)
            return (2);
        return (1);
    }
}


// Функция, которая проверяет, пересекаются ли два кубических сплайна, и если да, находит точки пересечения
int find_intersection(Coords *intersect_points, CubicSpline *spline1, CubicSpline *spline2, double *x_1, double *x_2,
                      int n, int m) {
    int curr = 0;
    double *roots = (double *) malloc(3 * sizeof(double));
    int roots_count;

    // Решаем уравнения на каждом интервале
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < m - 1; j++) {
            double left_border = max(x_1[i], x_2[j]);
            double right_border = min(x_1[i + 1], x_2[j + 1]);
            if (left_border <= right_border) {
                // ax^3 + bx^2 + cx + d
                // Находим разности коэффициентов для удобного решения по Кардано
                double a = spline1[i].d - spline2[j].d;
                double b = (spline1[i].c - 3 * spline1[i].d * x_1[i]);
                b -= (spline2[j].c - 3 * spline2[j].d * x_2[j]);

                double c = spline1[i].b - 2 * spline1[i].c * x_1[i] +
                           3 * spline1[i].d * pow(x_1[i], 2);
                c -= spline2[j].b - 2 * spline2[j].c * x_2[j] +
                     3 * spline2[j].d * pow(x_2[j], 2);

                double d = spline1[i].a - spline1[i].b * x_1[i] +
                           spline1[i].c * pow(x_1[i], 2) -
                           spline1[i].d * pow(x_1[i], 3);
                d -= spline2[j].a - spline2[j].b * x_2[j] +
                     spline2[j].c * pow(x_2[j], 2) -
                     spline2[j].d * pow(x_2[j], 3);

                if (a != 0) {
                    // кубический случай
                    // приводим уравнение для Кардано
                    b /= a;
                    c /= a;
                    d /= a;
                    roots_count = solve_cubic(roots, b, c, d);
                } else if (b != 0) {
                    // квадратичный случай
                    roots_count = solve_quadratic(roots, b, c, d);
                } else if (c != 0) {
                    // линейный случай
                    roots_count = 1;
                    roots[0] = -d / c;
                } else {
                    roots_count = 0;
                }

                for (int k = 0; k < roots_count; k++) {
                    if (roots[k] >= left_border && roots[k] <= right_border) {
                        // проверка совпадающих решений на концах областей определения функций
                        if (curr > 0 && intersect_points[curr - 1].x == roots[k])
                            continue;
                        // Записываем точки пересечения в массив
                        intersect_points[curr].x = roots[k];
                        intersect_points[curr].y = calculate_point(&spline1[i], roots[k], x_1[i]);
                        curr++;
                    }
                }
            }
        }
    }

    free(roots);

    // Общее количество точек
    return curr;
}

// Функция градиентного спуска нужна для нахождения минимального расстояния между сплайнами
void calculate_gradient(double *grad, CubicSpline *spline1_coef, CubicSpline *spline2_coef, double *x, double x10,
                        double x20) {
    // f(x1, x2) = (x1 - x2)^2 + (f1(x1) - f2(x2))^2
    // f(x1, x2) -> min
    // f' по x1 = 2(x1 - x2) + 2(f1(x1) - f2(x2)) * (b1 + 2c1(x1-x10) + 3d1(x1-x10)^2)
    // f' по x2 = -2(x1 - x2) + 2(f1(x1) - f2(x2)) * -(b2 + 2c2(x2-x20) + 3d2(x2-x20)^2)

    double f1 = calculate_point(spline1_coef, x[0], x10);
    double f2 = calculate_point(spline2_coef, x[1], x20);

    grad[0] = 2 * (x[0] - x[1]) + 2 * (f1 - f2) * (spline1_coef->b + 2 * spline1_coef->c * (x[0] - x10) +
                                                   3 * spline1_coef->d * pow(x[0] - x10, 2));

    grad[1] = -2 * (x[0] - x[1]) + 2 * (f1 - f2) * -(spline2_coef->b + 2 * spline2_coef->c * (x[1] - x20) +
                                                     3 * spline2_coef->d * pow(x[1] - x20, 2));
}


// Функция, которая находит минимальное расстояние между двумя кубическими сплайнами

double find_min_distance(CubicSpline *spline1, CubicSpline *spline2, double *x_1, double *x_2, int n, int m) {
    // f(x1, x2) = (x1 - x2) ^ 2 + (f1(x1) - f2(x2)) ^ 2
    // f(x1, x2) -> min
    // f' по x1 = 2(x1 - x2) + 2(f1(x1) - f2(x2)) * (b1 + 2c1(x1-x10) + 3d1(x1-x10)^2)
    // f' по x2 = -2(x1 - x2) + 2(f1(x1) - f2(x2)) * -(b2 + 2c2(x2-x20) + 3d2(x2-x20)^2)

    double grad[2];
    double x_prev[2], x[2];
    double lambda = 0.1;
    double minimum = INFINITY;
    short flag = 1;
    int counts = 0;

    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < m - 1; j++) {
            // Решение
            x_prev[0] = x_1[i];
            x_prev[1] = x_2[j];

            while (1) {
                calculate_gradient(grad, &spline1[i], &spline2[j], x_prev, x_1[i],
                                   x_2[j]);
                x[0] = x_prev[0] - lambda * grad[0];
                x[1] = x_prev[1] - lambda * grad[1];

                if ((fabs(x[0] - x_prev[0]) < EPSILON && fabs(x[1] - x_prev[1]) < EPSILON))
                    break;

                if (x[0] < x_1[i] || x[0] > x_1[i + 1] ||
                    x[1] < x_2[j] || x[1] > x_2[j + 1])
                    break;

                x_prev[0] = x[0];
                x_prev[1] = x[1];
            }

            if (x[0] >= x_1[i] && x[0] <= x_1[i + 1] &&
                x[1] >= x_2[j] && x[1] <= x_2[j + 1]) {
                // Находим расстояние между точками
                double res = sqrt(pow(x[0] - x[1], 2) + pow(calculate_point(&spline1[i], x[0], x_1[i]) -
                                                       calculate_point(&spline2[j], x[1], x_2[j]), 2));
                if (res < EPSILON) {
                    flag = 0;
                    counts++;
                }
                minimum = min(minimum, res);
            }
        }
    }

    return flag ? -1 : minimum;
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


void sort(double *arr_x, double *arr_y, int n, int k) {
    int noSwap;
    double tmp;
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

    // Если срочно нужно посчитать сплайны
//    int n = 4;
//    // инициализация координат x и y точек сплайнов
//    double x1[] = {1, 2, 4, 7};
//    double y1[] = {2, 3, 1, 4};
//
//    int m = 4;
//    double x2[] = {2, 4, 7, 11};
//    double y2[] = {3, 4, 2, 5};
//    printf("------------------------------------------------------------------------------------------\n");

    // сортировка сплайнов
    sort(x1, y1, n, 1);
    sort(x2, y2, m, 2);

    // создание сплайнов
    CubicSpline spline10[n];
    CubicSpline spline20[m];

    // вычисление коэффициентов
    compute_spline_coefficients(x1, y1, n, spline10);
    compute_spline_coefficients(x2, y2, m, spline20);

    //Здесь мы избавляемся от тех самых косячных последних элементов
    CubicSpline spline1[n - 1];
    CubicSpline spline2[m - 1];

    for (int i = 0; i < n - 1; i++) {
        spline1[i] = spline10[i];
    }

    for (int i = 0; i < m - 1; i++) {
        spline2[i] = spline20[i];
    }
    printf("Coefficients of spline_1 and spline_2:\n");
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

    //Массив для точек пересечения сплайнов
    Coords *intersect_points = (Coords *) malloc((n - 1) * (m - 1) * 3 * sizeof(Coords));
    // проверка на пересечение сплайнов
    // ВОТ НАД ЭТИМ НАДО РАБОТАТЬ
    int cross_p = find_intersection(intersect_points, spline1, spline2, x1, x2, n, m);
    if (cross_p == 0) {
        printf("Spline_1 and spline_2 don't intersect\n");
        double min_dist = find_min_distance(spline1, spline2, x1, x2, n, m);
        if (min_dist != -1) {
            printf("The minimum distance between spline 1 and spline 2 is %f.\n", min_dist);
        } else {
            printf("The minimum distance can't be found.");
        }
    } else {
        printf("Intersection points are: \n");
        for (int i = 0; i < cross_p; i++){
            printf("x = %lf; y = %lf\n", intersect_points[i].x, intersect_points[i].y);
        }
    }

    free(intersect_points);
    return 0;
}