#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define MAX_POINTS 1000

/* Кубический сплайн для интерполяции */
typedef struct {
    double a, b, c, d, x;
} Spline;

/* Вычисление коэффициентов сплайна */
void compute_spline_coefficients(double *x, double *y, int n, Spline *splines) {
    double *alpha = (double *) malloc(n * sizeof(double));
    double *beta = (double *) malloc(n * sizeof(double));
    double *gamma = (double *) malloc(n * sizeof(double));
    double *delta = (double *) malloc(n * sizeof(double));

    for (int i = 1; i < n; i++) {
        double dx1 = x[i] - x[i - 1];
        double dy1 = y[i] - y[i - 1];
        double dx2 = x[i + 1] - x[i];
        double dy2 = y[i + 1] - y[i];

        // Проверяем, что разность x-координат двух точек не равна нулю
        if (dx1 == 0 || dx2 == 0) {
            // Не можем вычислить коэффициенты сплайна, устанавливаем их в ноль
            splines[i].a = 0;
            splines[i].b = 0;
            splines[i].c = 0;
            splines[i].d = 0;
            splines[i].x = x[i];
            continue;
        }

        double hi = dx1;
        double hi1 = dx2;
        double dydx1 = dy1 / dx1;
        double dydx2 = dy2 / dx2;
        double diff_dydx = dydx2 - dydx1;
        double bi = 2 * (hi + hi1);
        double ci = hi1;
        double di = 6 * (diff_dydx / (hi1 + hi));

        alpha[i] = hi1 / bi;
        beta[i] = -(hi + hi1) / bi;
        gamma[i] = 1;
        delta[i] = di / bi;
    }

    // Прямой ход метода прогонки
    double *c = (double *) malloc(n * sizeof(double));
    double *l = (double *) malloc(n * sizeof(double));
    double *mu = (double *) malloc(n * sizeof(double));
    double *z = (double *) malloc(n * sizeof(double));

    l[0] = 1;
    mu[0] = 0;
    z[0] =
    c[0] = 0;

    for (int i = 1; i < n; i++) {
        l[i] = gamma[i] - alpha[i] * mu[i - 1];
        mu[i] = c[i] / l[i];
        z[i] = (delta[i] - alpha[i] * z[i - 1]) / l[i];
    }

// Обратный ход метода прогонки
    c[n - 1] = z[n - 1];

    for (int i = n - 2; i >= 0; i--) {
        c[i] = z[i] - mu[i] * c[i + 1];
    }

// Вычисление коэффициентов сплайнов
    for (int i = 1; i < n; i++) {
        double dx = x[i] - x[i - 1];
        splines[i].a = (c[i] - c[i - 1]) / (6 * dx);
        splines[i].b = 0.5 * c[i];
        splines[i].c = (y[i] - y[i - 1]) / dx - (2 * dx * c[i] + dx * c[i - 1]) / 6;
        splines[i].d = y[i - 1];
        splines[i].x = x[i];
    }

    free(alpha);
    free(beta);
    free(gamma);
    free(delta);
    free(c);
    free(l);
    free(mu);
    free(z);
}

/* Поиск минимального расстояния между сплайнами */
double minimum_distance_between_splines(Spline *s1, Spline *s2) {
// Если сплайны совпадают, расстояние равно нулю
    if (s1->x == s2->x && s1->a == s2->a && s1->b == s2->b && s1->c == s2->c && s1->d == s2->d) {
        return 0.0;
    }
    // Находим минимальное расстояние между сплайнами, решая квадратное уравнение
    double a = s1->a - s2->a;
    double b = s1->b - s2->b;
    double c = s1->c - s2->c;
    double d = s1->d - s2->d;

    double A = 3 * a * c - b * b;
    double B = 2 * b * b - 4 * a * d + 2 * a * s2->d - 2 * b * s2->c;
    double C = a * d - b * s2->d + c * s2->c;

    double delta = B * B - 4 * A * C;

    if (delta < 0) {
        // Нет реальных корней, минимальное расстояние находится на концах отрезков
        double dist1 = fabs(s1->d - s2->d);
        double dist2 = fabs((s1 + 1)->d - s2->d);
        double dist3 = fabs(s1->d - (s2 + 1)->d);
        double dist4 = fabs((s1 + 1)->d - (s2 + 1)->d);
        return fmin(fmin(dist1, dist2), fmin(dist3, dist4));
    } else {
        // Решаем квадратное уравнение
        double t1 = (-B + sqrt(delta)) / (2 * A);
        double t2 = (-B - sqrt(delta)) / (2 * A);

        // Находим точки на сплайнах, на которых достигается минимальное расстояние
        double x1 = s1->x + t1;
        double x2 = s1->x + t2;

        // Проверяем, находятся ли точки на отрезках сплайнов
        if (x1 < s1->x || x1 > (s1 + 1)->x || x1 < s2->x || x1 > (s2 + 1)->x) {
            x1 = NAN;
        }

        if (x2 < s1->x || x2 > (s1 + 1)->x || x2 < s2->x || x2 > (s2 + 1)->x) {
            x2 = NAN;
        }

        // Вычисляем расстояния до точек
        double dist1 = (x1 != NAN) ? fabs(spline_value(s1, x1) - spline_value(s2, x1)) : INFINITY;
        double dist2 = (x2 != NAN) ? fabs(spline_value(s1, x2) - spline_value(s2, x2)) : INFINITY;

        // Возвращаем минимальное расстояние
        return fmin(dist1, dist2);
    }
}

/* Проверка на пересечение сплайнов */
bool splines_intersect(Spline *s1, Spline *s2) {
// Если сплайны совпадают, они пересекаются
    if (s1->x == s2->x && s1->a == s2->a && s1->b == s2->b && s1->c == s2->c && s1->d == s2->d) {
        return true;
    }

// Находим точки пересечения сплайнов, решая систему уравнений
    double a1 = s1->a;
    double b1 = s1->b;
    double c1 = s1->c - s2->c;
    double d1 = s1->d - s2->d;

    double a2 = s2->a;
    double b2 = s2->b;
    double c2 = 0;
    double d2 = 0;

    double t1, t2;

    if (fabs(a1) < 1e-10) {
        t1 = -d1 / c1;
        t2 = (d2 - b2 * t1) / (a2 * t1 + c2);
    } else if (fabs(a2) < 1e-10) {
        t2 = -d2 / c2;
        t1 = (d1 - b1 * t2) / (a1 * t2 + c1);
    } else {
        double k = a2 / a1;
        double q = (c2 - k * c1) / (b2 - k * b1);
        double p = (d2 - k * d1) / (b2 - k * b1);

        double A = k * k * a1 - a2;
        double B = 2 * k * q * a1 - b2 + k * b1;
        double C = q * q * a1 - c2 + k * c1;
        double D = p * p * a1 - d2 + k * d1;

        double t1_1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
        double t1_2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
        double t2_1 = q - k * t1_1;
        double t2_2 = q - k * t1_2;

        t1 = isnan(t1_1) ? t1_2 : (isnan(t1_2) ? t1_1 : fmin(t1_1, t1_2));
        t2 = isnan(t2_1) ? t2_2 : (isnan(t2_2) ? t2_1 : fmin(t2_1, t2_2));
    }

// Проверяем, находятся ли точки пересечения на отрезках сплайнов
    if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
        return true;
    } else {
        return false;
    }
}

/* Функция для вычисления значения сплайна в точке */
int spline_value(Spline *s, double x) {
    double dx = x - s->x;
    return s->a + s->b * dx + s->c * dx * dx + s->d * dx * dx * dx;
}

/* Функция для нахождения минимального расстояния между сплайнами */
double find_min_distance(Spline *s1, Spline *s2) {
    double x, y1, y2, dist;
    double min_dist = INFINITY;

// Ищем минимальное расстояние между точками на сплайнах
    for (x = s1->x; x <= s1->x + s1->h; x += 0.01) {
        y1 = spline_value(s1, x);
        y2 = spline_value(s2, x);
        dist = sqrt(pow(x - s1->x, 2) + pow(y1 - y2, 2));
        if (dist < min_dist) {
            min_dist = dist;
        }
    }

// Ищем пересечения сплайнов
    if (check_intersection(s1, s2)) {
        min_dist = 0;
    }

    return min_dist;
}

void hello() {
    double *x, *y;
    Spline *splines;
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
    while (n0 <= 1 || n0 != (int) n0) {
        printf("Please, enter a positive integer starting from 2:\n");
        scanf("%lf", &n0);
    }
    int n = (int) n0;

//   Вводим x
    printf("Enter the function arguments(x):\n");
    // Выделяем память для массивов координат
    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        scanf("%lf", &x[i]);
    }

//  Вводим y
    printf("Enter the function values(y):\n");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &y[i]);
    }
    // Создаем сплайны
    splines = cubic_spline(x, y, n);

// Ищем минимальное расстояние между сплайнами
    double min_dist = INFINITY;
    for (int i = 0; i < n - 1; i++) {
        int j;
        for (j = i + 1; j < n; j++) {
            double dist = find_min_distance(&splines[i], &splines[j]);
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
    }

// Выводим результат
    if (min_dist == INFINITY) {
        printf("No splines found.\n");
    } else {
        printf("Minimum distance between splines: %.2f\n", min_dist);
    }

// Освобождаем память
    free(x);
    free(y);
    free(splines);
}


int main() {
    hello();
    return 0;
}
