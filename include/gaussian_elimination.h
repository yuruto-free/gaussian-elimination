#ifndef GAUSSIAN_ELIMINATION_H__
#define GAUSSIAN_ELIMINATION_H__

#include <stdint.h>

#define GE_OK (0)
#define GE_NG (1)

/**
 * \brief ガウスの消去法
 * \param[in]    dim    次元数
 * \param[in]    matrix 係数行列
 *                      サイズ：dim * dim
 * \param[inout] vec    右辺ベクトル
 *                      サイズ：dim
 * \retval       GE_OK  正常終了
 * \retval       GE_NG  異常終了
 */
int32_t GE_gauss_solver(int32_t dim, double *matrix, double *vec);

#endif
