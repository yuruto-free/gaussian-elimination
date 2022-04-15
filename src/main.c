#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include "gaussian_elimination.h"
#define MAX_TEST_DATA (3)

struct data_t {
    int32_t dim;
    double *matrix;
    double *vec;
    double *ans;
};

double matrix1[9] = {
    2.0, 1.0, -2.0, 
    1.0, 1.0, -1.0, 
    1.0, -2.0, 3.0,
};
double vec1[3] = {
    1.0, 4.0, -1.0,
};
double ans1[3] = {
    1.0, 7.0, 4.0,
};

double matrix2[9] = {
    0.0, 0.0, 2.0, 
    3.0, 3.0, 1.0, 
    -2.0, 2.0, 0.0,
};
double vec2[3] = {
    1.0, 2.0, 3.0,
};
double ans2[3] = {
    -0.5, 1.0, 0.5, 
};

double matrix3[16] = {
    3.0, -3.0, 4.0, 2.0,
    7.0, -5.0, -2.0, -1.0,
    3.0, 7.0, 9.0, 1.0,
    1.0, 2.0, 3.0, 4.0,
};
double vec3[4] = {
    17.0, -13.0, 48.0, 30.0,
};
double ans3[4] = {
    1.0, 2.0, 3.0, 4.0,
};

struct data_t TEST_DATA[MAX_TEST_DATA] = {
    {3, matrix1, vec1, ans1},
    {3, matrix2, vec2, ans2},
    {4, matrix3, vec3, ans3},
};

static double calculate_error(const struct data_t *target) {
    int32_t i;
    int32_t dim;
    double sum;
    dim = target->dim;
    sum = 0.0;

    for (i = 0; i < dim; i++) {
        sum += fabs(target->vec[i] - target->ans[i]);
    }

    return sum;
}

static void print_vec(const struct data_t *target) {
    int32_t i;
    int32_t dim;
    dim = target->dim;

    for (i = 0; i < dim; i++) {
        printf(" %.5f", target->vec[i]);
    }
    printf("\n");
}

int main(int argc, char **argv) {
    struct data_t target;
    double err;
    int32_t ret;
    int i;

    for (i = 0; i < (int)MAX_TEST_DATA; i++) {
        memcpy(&target, &TEST_DATA[i], sizeof(struct data_t));
        printf("Dim: %d\n", target.dim);
        ret = GE_gauss_solver(target.dim, target.matrix, target.vec);
        if ((int32_t)GE_OK != ret) {
            goto EXIT_MAIN;
        }
        err = calculate_error((const struct data_t *)&target);
        print_vec((const struct data_t *)&target);
        printf("error: %.13e\n", err);
        printf("\n");
    }

EXIT_MAIN:

    return 0;
}