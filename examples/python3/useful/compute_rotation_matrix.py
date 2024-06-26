import numpy as np

def rotation_matrix_x(angle):
    rotation = np.array([[1, 0, 0],
                         [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle), np.cos(angle)]])
    return rotation


def rotation_matrix_y(angle):
    rotation = np.array([[np.cos(angle), 0, np.sin(angle)],
                         [0, 1, 0],
                         [-np.sin(angle), 0, np.cos(angle)]])
    return rotation


def rotation_matrix_z(angle):
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0],
                         [0, 0, 1]])
    return rotation

def compute_rotation_matrix(x, y, z):
    rotation = np.dot(rotation_matrix_z(z), np.dot(rotation_matrix_y(y), rotation_matrix_x(x)))
    return rotation
