#include "gaussian_elimination.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#define RETURN_OK_GE (0)
#define RETURN_NG_GE (1)
#define MEPS_GE (1e-10)

/*
 * \brief ピボット選択
 * \param[in]  dim    次元数
 * \param[in]  col    対象列
 * \param[in]  matrix 係数行列
 * \param[out] pivot  ピボット値
 * \retval RETURN_OK_GE 正常終了
 * \retval RETURN_NG_GE 異常終了
 */
static int32_t pivot_selection(int32_t dim, int32_t col, double *matrix,
                               int32_t *pivot);
/*
 * \brief 入れ替え処理
 * \param[in]    dim    次元数
 * \param[in]    row    入れ替え元の行
 * \param[in]    pivot  入れ替え先の行
 * \param[inout] matrix 係数行列
 * \param[inout] vec    右辺ベクトル
 * \retval RETURN_OK_GE 正常終了
 * \retval RETURN_NG_GE 異常終了
 */
static int32_t swap(int32_t dim, int32_t row, int32_t pivot, double *matrix,
                    double *vec);

int32_t GE_gauss_solver(int32_t dim, double *matrix, double *vec) {
    int32_t ret = (int32_t)GE_NG;
    int32_t func_val;
    int32_t pivot, idx;
    int32_t row, col;
    double diag, scale, sum;

    if ((NULL != matrix) && (NULL != vec)) {
        /* 前進消去 */
        for (idx = 0; idx < dim - 1; idx++) {
            /* ピボット選択 */
            func_val = pivot_selection(dim, idx, matrix, &pivot);
            if ((int32_t)RETURN_OK_GE != func_val) {
                goto EXIT_GE_SOLVER;
            }
            /* 対角成分の絶対値が最大でない場合、行を入れ替える */
            if (idx != pivot) {
                func_val = swap(dim, idx, pivot, matrix, vec);
                if ((int32_t)RETURN_OK_GE != func_val) {
                    goto EXIT_GE_SOLVER;
                }
            }
            /* 対角成分の最大値を取得 */
            diag = matrix[idx * dim + idx];

            if (fabs(diag) < (double)MEPS_GE) {
                goto EXIT_GE_SOLVER;
            }
            for (row = idx + 1; row < dim; row++) {
                scale = matrix[row * dim + idx] / diag;

                for (col = idx; col < dim; col++) {
                    matrix[row * dim + col] -= matrix[idx * dim + col] * scale;
                }
                vec[row] -= vec[idx] * scale;
            }
        }

        /* 後退代入 */
        for (row = dim - 1; row >= 0; row--) {
            sum = vec[row];

            for (col = row + 1; col < dim; col++) {
                sum -= matrix[row * dim + col] * vec[col];
            }
            vec[row] = sum / matrix[row * dim + row];
        }
        ret = (int32_t)GE_OK;
    }
EXIT_GE_SOLVER:

    return ret;
}

static int32_t pivot_selection(int32_t dim, int32_t col, double *matrix,
                               int32_t *pivot) {
    int32_t ret = (int32_t)RETURN_NG_GE;
    int32_t i;
    int32_t p;
    double max_val, val;

    if ((NULL != matrix) && (NULL != pivot)) {
        /* 対角成分の絶対値が最大値となると仮定 */
        max_val = fabs(matrix[col * dim + col]);
        p       = col;

        /* 対角成分より下側の値のうち、絶対が最大となるものを探索 */
        for (i = col + 1; i < dim; i++) {
            val = fabs(matrix[i * dim + col]);

            if (max_val < val) {
                max_val = val;
                p       = i;
            }
        }
        (*pivot) = p;
        ret      = (int32_t)RETURN_OK_GE;
    }

    return ret;
}

static int32_t swap(int32_t dim, int32_t row, int32_t pivot, double *matrix,
                    double *vec) {
    int32_t ret = (int32_t)RETURN_NG_GE;
    int32_t col;
    double tmp;

    if ((NULL != matrix) && (NULL != vec)) {
        /* 行の入れ替え */
        for (col = 0; col < dim; col++) {
            tmp                       = matrix[row * dim + col];
            matrix[row * dim + col]   = matrix[pivot * dim + col];
            matrix[pivot * dim + col] = tmp;
        }
        /* 右辺ベクトルの入れ替え */
        tmp        = vec[row];
        vec[row]   = vec[pivot];
        vec[pivot] = tmp;
        ret        = (int32_t)RETURN_OK_GE;
    }

    return ret;
}