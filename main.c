#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Константа, определяющая точность вычислений
#define EPSILON 1e-6

// Структура, описывающая кубический сплайн
typedef struct {
    double a, b, c, d; // Коэффициенты для каждого интервала
    double x; // Левый конец интервала
} CubicSpline;

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

// Функция, которая проверяет, пересекаются ли два кубических сплайна, и если да, находит точки пересечения
int find_intersection(CubicSpline *spline1, CubicSpline *spline2, double *x_intersect, double *y_intersect, int n) {
    double t1, t2, t3, t4, a, b, c, delta, t;
    int i, intersect_count = 0;
    // Проверяем пересечение каждого интервала
    for (i = 0; i < n - 1; i++) {
        t1 = spline1[i].a - spline2[i].a + spline2[i].b * spline1[i].x - spline1[i].b * spline2[i].x;
        t2 = spline1[i].c - spline2[i].c;
        t3 = spline1[i].d - spline2[i].d;
        t4 = spline1[i + 1].a - spline2[i + 1].a + spline2[i + 1].b * spline1[i + 1].x -
             spline1[i + 1].b * spline2[i + 1].x;
        a = -t1 + t4;
        b = 3 * t1 - 2 * t2 + t4;
        c = -2 * t1 + t2 + t3 - t4;
        delta = b * b - 4 * a * c;
        if (delta >= 0) {
            // Если уравнение имеет корни, проверяем их значения
            if (fabs(a) < EPSILON) {
                t = -c / b;
                if (t >= 0 && t <= 1) {
                    x_intersect[intersect_count] = spline1[i].x + t * (spline1[i + 1].x - spline1[i].x);
                    y_intersect[intersect_count] =
                            spline1[i].a + t * (spline1[i].b + t * (spline1[i].c + t * spline1[i].d));
                    intersect_count++;
                }
            } else {
                t = (-b + sqrt(delta)) / (2 * a);
                if (t >= 0 && t <= 1) {
                    x_intersect[intersect_count] = spline1[i].x + t * (spline1[i + 1].x - spline1[i].x);
                    y_intersect[intersect_count] =
                            spline1[i].a + t * (spline1[i].b + t * (spline1[i].c + t * spline1[i].d));
                    intersect_count++;
                }
                t = (-b - sqrt(delta)) / (2 * a);
                if (t >= 0 && t <= 1) {
                    x_intersect[intersect_count] = spline1[i].x + t * (spline1[i + 1].x - spline1[i].x);
                    y_intersect[intersect_count] =
                            spline1[i].a + t * (spline1[i].b + t * (spline1[i].c + t * spline1[i].d));
                    intersect_count++;
                }
            }
        }
    }
    return intersect_count;
}


// Функция, которая находит минимальное расстояние между двумя кубическими сплайнами
double find_min_distance(CubicSpline *spline1, CubicSpline *spline2, int n, int m) {
    double min_distance = INFINITY;
    int i, j;
    // Вычисляем расстояние между каждой парой точек на обоих сплайнах
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < m - 1; j++) {
            double x1 = spline1[i].x;
            double y1 = spline1[i].a + spline1[i].b * (x1 - spline1[i].x) + spline1[i].c * pow(x1 - spline1[i].x, 2) +
                        spline1[i].d * pow(x1 - spline1[i].x, 3);
            double x2 = spline2[j].x;
            double y2 = spline2[j].a + spline2[j].b * (x2 - spline2[j].x) + spline2[j].c * pow(x2 - spline2[j].x, 2) +
                        spline2[j].d * pow(x2 - spline2[j].x, 3);
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
    double y = spline[i].a + spline[i].b * h + spline[i].c * pow(h, 2)+spline[i].d * pow(h, 3);
    return y;
}


void sort(double *arr_x, double *arr_y, int n){
    int noSwap;
    int tmp;
    for (int i = n - 1; i >= 0; i--)
    {
        noSwap = 1;
        for (int j = 0; j < i; j++)
        {
            if (arr_x[j] > arr_x[j + 1])
            {
                tmp = arr_x[j];
                arr_x[j] = arr_x[j + 1];
                arr_x[j + 1] = tmp;
                noSwap = 0;
                tmp= arr_y[j];
                arr_y[j]=arr_y[j+1];
                arr_y[j+1]=tmp;
            }
        }
        if (noSwap == 1)
            break;
    }
    printf("Sorted value of x and y\n");
    for (int i=0; i<n;i++){
        printf("%lf %lf\n",arr_x[i],arr_y[i]);
    }

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
           "Let's enter the number of our points to splain 1:\n", name);

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

    // Также делаем и для второго сплайна
    printf("Let's enter the number of our points to splain 2:\n");

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

    // создание сплайнов
    sort(x1,y1,n);
    sort(x2,y2,m);
    CubicSpline spline1[n];
    CubicSpline spline2[m];
    compute_spline_coefficients(x1, y1, n, spline1);
    compute_spline_coefficients(x2, y2, n, spline2);

    //Выводим коэффициенты
    for (int i = 0; i < n; i++) {
        printf("Spline%d: a = %lf, b = %lf, c = %lf, d = %lf, x = %lf\n", i, spline1[i].a, spline1[i].b, spline1[i].c,
               spline1[i].d, spline1[i].x);
    }
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("Spline%d: a = %lf, b = %lf, c = %lf, d = %lf, x = %lf\n", i, spline2[i].a, spline2[i].b, spline2[i].c,
               spline2[i].d, spline2[i].x);
    }

    // Спрашиваем у пользователя, в каком сплайне мы хотим вычислить значение точки сплайна
    // И вычисляем
    double new_x;
    int spline_num;
    double new_y;
    printf("Enter the point which function value you want to know\n");
    scanf("%lf",&new_x);
    printf("Enter the spline number which you would like to work on\n");
    scanf("%d",&spline_num);
    while (spline_num!=1 && spline_num!=2){
        printf("The spline number can be only 1 or 2\n");
        scanf("%d",&spline_num);
    }
    if(spline_num == 1){
        new_y=evaluate_cubic_spline( new_x, x1, spline1, n);
    }
    else {
        new_y=evaluate_cubic_spline( new_x, x2, spline2, m);
    }

    printf("The value of function in point %lf = %lf\n", new_x, new_y);

    // проверка совпадения сплайнов
    int is_same_spline = is_cubic_spline_same(spline1, spline2, n, m);
    if (is_same_spline) {
        printf("Spline 1 and spline 2 are the same.\n");
        return 0;
    }

    // проверка на пересечение сплайнов
    double x_intersect[2], y_intersect[2];
    int intersect_count = find_intersection(spline1, spline2, x_intersect, y_intersect, n);

    if (intersect_count > 0) {
        // Если сплайны пересекаются, выводим координаты точек пересечения
        printf("Spline_1 and spline_2 intersect at:\n");
        for (int i = 0; i < intersect_count; i++) {
            printf("(%f, %f)\n", x_intersect[i], y_intersect[i]);
        }
    } else {
        // Если сплайны не пересекаются, находим минимальное расстояние между ними
        double min_distance = find_min_distance(spline1, spline2, n, m);
        printf("The minimum distance between spline 1 and spline 2 is %f.\n", min_distance);
    }

    return 0;
}
