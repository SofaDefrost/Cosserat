import numpy as np
from compute_rotation_matrix import *


def piecewise_logmap1(curv_abs, g_x):
    # xi_hat = np.zeros((4, 4), dtype=float)

    theta = compute_theta(curv_abs, g_x)

    if theta == 0.0:
        xi_hat = 1.0 / curv_abs * (g_x - np.identity(4))
    else:
        t0 = curv_abs * theta
        t1 = np.sin(t0)
        t2 = np.cos(t0)
        t3 = 2 * t1 * t2
        t4 = 1 - 2 * t1 ** 2
        t5 = t0 * t4

        gp2 = np.dot(g_x, g_x)
        gp3 = np.dot(gp2, g_x)

        xi_hat = 1.0 / curv_abs * (0.125 * (1.0 / np.sin(t0 / 2.0) ** 3) * np.cos(t0 / 2.0) *
                                   ((t5 - t1) * np.identity(4) - (t0 * t2 + 2 * t5 - t1 - t3) * g_x +
                                    (2 * t0 * t2 + t5 - t1 - t3) * gp2 - (t0 * t2 - t1) * gp3))

    xci = np.array([xi_hat[2, 1], xi_hat[0, 2], xi_hat[1, 0], xi_hat[0, 3], xi_hat[1, 3], xi_hat[2, 3]])

    return xci


def piecewise_logmap2(curv_abs, gX):
    theta = compute_theta(curv_abs, gX)
    I4 = np.identity(4)
    xi_hat = np.zeros((4, 4))

    csc_theta = 1.0 / np.sin(curv_abs * theta / 2.0)
    sec_theta = 1.0 / np.cos(curv_abs * theta / 2.0)
    cst = (1.0 / 8) * (csc_theta ** 3) * sec_theta
    x_theta = curv_abs * theta
    cos_2Xtheta = np.cos(2.0 * x_theta)
    cos_Xtheta = np.cos(x_theta)
    sin_2Xtheta = np.sin(2.0 * x_theta)
    sin_Xtheta = np.sin(x_theta)

    if theta <= np.finfo(float).eps:
        xi_hat = I4
    else:
        xi_hat = cst * ((x_theta * cos_2Xtheta - sin_Xtheta) * I4 -
                        (x_theta * cos_Xtheta + 2.0 * x_theta * cos_2Xtheta - sin_Xtheta - sin_2Xtheta) * gX +
                        (2.0 * x_theta * cos_Xtheta + x_theta * cos_2Xtheta - sin_Xtheta - sin_2Xtheta) * (
                            np.dot(gX, gX)) -
                        (x_theta * cos_Xtheta - sin_Xtheta) * (np.dot(np.dot(gX, gX), gX)))

    xi = np.array([xi_hat[2, 1], xi_hat[0, 2], xi_hat[1, 0], xi_hat[0, 3], xi_hat[1, 3], xi_hat[2, 3]])
    return xi


def compute_theta(x, gX):
    Tr_gx = np.trace(gX)
    if x <= np.finfo(float).eps:
        theta = 0.0
    else:
        theta = (1.0 / x) * np.arccos((Tr_gx / 2.0) - 1)

    return theta


if __name__ == '__main__':
    _curv_abs = 1
    angle_y = 20 * np.pi / 180


    R = rotation_matrix_y(angle_y)

    print(f'The rotation matrix is: \n {R}')

    _g_x = np.zeros((4, 4), dtype=float)
    for i in range(3):
        for j in range(3):
            _g_x[i][j] = R[i][j]

    _g_x[0][3] = _curv_abs  # to deploy the beam node and the rest part of transform is equal to null
    _g_x[3][3] = 1  # The homogeneous matrix

    xci = piecewise_logmap1(_curv_abs, _g_x)

    print(f'The piecewise logmap is: {xci}')

    # 0.29241528
