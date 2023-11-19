import math


def looptest(int, loopnum):
    vertex_loop = []
    loops = []
    spring_set = []
    spring = []
    Ks = 10e12
    Kd = 0
    indices1 = []
    indices2 = []
    lengths = []

    for loop_id in range(0, loopnum):
        for v_id in range(0, 45):
            if v_id == 0:
                vertex_loop.append(
                    [int * loop_id + 1, 7.5 * math.cos(v_id / 25 * math.pi), 7.5 * math.sin(v_id / 25 * math.pi)])
            if 0 < v_id < 25:
                vertex_loop.append(
                    [int * loop_id + 1, 7.5 * math.cos(v_id / 25 * math.pi), 7.5 * math.sin(v_id / 25 * math.pi)])
                L = (math.dist(vertex_loop[v_id], vertex_loop[v_id - 1]))
                spring.append([v_id - 1, v_id, Ks, Kd, L])
                indices1.append(v_id - 1)
                indices2.append(v_id)
                lengths.append(L)
            elif 25 <= v_id < 28:
                vertex_loop.append([int * (loop_id) + 1, 7.5 * math.cos(math.pi), 7.5 * math.sin(math.pi) - v_id + 25])
                L = (math.dist(vertex_loop[v_id], vertex_loop[v_id - 1]))
                spring.append([v_id - 1, v_id, Ks, Kd, L])
                indices1.append(v_id - 1)
                indices2.append(v_id)
                lengths.append(L)
            elif v_id >= 28 and v_id < 42:
                vertex_loop.append([int * (loop_id) + 1, 7.5 * math.cos(math.pi) + 1 + v_id - 28, -4])
                L = (math.dist(vertex_loop[v_id], vertex_loop[v_id - 1]))
                spring.append([v_id - 1, v_id, Ks, Kd, L])
                indices1.append(v_id - 1)
                indices2.append(v_id)
                lengths.append(L)
            elif v_id >= 42:
                vertex_loop.append([int * (loop_id) + 1, -7.5 + 15, -2])
                L = (math.dist(vertex_loop[v_id], vertex_loop[v_id - 1]))
                spring.append([v_id - 1, v_id, Ks, Kd, L])
                indices1.append(v_id - 1)
                indices2.append(v_id)
                lengths.append(L)
            if v_id == 43:
                L = (math.dist(vertex_loop[v_id], vertex_loop[0]))
                spring.append([v_id, 0, Ks, Kd, L])
                indices1.append(v_id - 1)
                indices2.append(0)
                lengths.append(L)
        spring_set.append(spring)
        loops.append(vertex_loop)
        vertex_loop = []
        spring = []
    return loops, spring_set, indices1, indices2, lengths

    # fiber.addObject('MeshTopology', src='@loader', name='container')
    # fiber.addObject('MechanicalObject', name='lines', template='Vec3', showObject=True, showObjectScale=1)


# def generate_loop_data(loop_num, loop_spacing):
#     loops = []
#     spring_set = []
#     indices1 = []
#     indices2 = []
#     lengths = []
#     Ks = 10e12
#     Kd = 0
#
#     def generate_vertex(loop_id, v_id):
#         if v_id == 0:
#             return [loop_id * loop_spacing + 1, 7.5 * math.cos(v_id / 25 * math.pi), 7.5 * math.sin(v_id / 25 * math.pi)]
#         elif 0 < v_id < 25:
#             L = math.dist(loops[loop_id][v_id], loops[loop_id][v_id - 1])
#             spring_set[loop_id].append([v_id - 1, v_id, Ks, Kd, L])
#             indices1.append(v_id - 1)
#             indices2.append(v_id)
#             lengths.append(L)
#             return [loop_id * loop_spacing + 1, 7.5 * math.cos(v_id / 25 * math.pi), 7.5 * math.sin(v_id / 25 * math.pi)]
#         elif 25 <= v_id < 28:
#             L = math.dist(loops[loop_id][v_id], loops[loop_id][v_id - 1])
#             spring_set[loop_id].append([v_id - 1, v_id, Ks, Kd, L])
#             indices1.append(v_id - 1)
#             indices2.append(v_id)
#             lengths.append(L)
#             return [loop_id * loop_spacing + 1, 7.5 * math.cos(math.pi), 7.5 * math.sin(math.pi) - v_id + 25]
#         elif 28 <= v_id < 42:
#             L = math.dist(loops[loop_id][v_id], loops[loop_id][v_id - 1])
#             spring_set[loop_id].append([v_id - 1, v_id, Ks, Kd, L])
#             indices1.append(v_id - 1)
#             indices2.append(v_id)
#             lengths.append(L)
#             return [loop_id * loop_spacing + 1, 7.5 * math.cos(math.pi) + 1 + v_id - 28, -4]
#         else:
#             L = math.dist(loops[loop_id][v_id], loops[loop_id][v_id - 1])
#             spring_set[loop_id].append([v_id - 1, v_id, Ks, Kd, L])
#             indices1.append(v_id - 1)
#             indices2.append(v_id)
#             lengths.append(L)
#             return [loop_id * loop_spacing + 1, -7.5 + 15, -2]
#
#     for loop_id in range(loop_num):
#         vertex_loop = [generate_vertex(loop_id, v_id) for v_id in range(45)]
#         L = math.dist(vertex_loop[43], vertex_loop[0])
#         spring_set[loop_id].append([43, 0, Ks, Kd, L])
#         indices1.append(43)
#         indices2.append(0)
#         lengths.append(L)
#         loops.append(vertex_loop)
#
#     return loops, spring_set, indices1, indices2, lengths
#
# loop_num = 5
# loop_spacing = 10
# loops_data, spring_data, indices1, indices2, lengths = generate_loop_data(loop_num, loop_spacing)
#
# print(loops_data)
# print(spring_data)
# print(indices1)
# print(indices2)
# print(lengths)
