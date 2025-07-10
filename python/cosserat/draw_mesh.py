# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Define the vertices of the tetrahedron
# vertices = [
#     [0, 0, 0],
#     [0.6, 0, 0],
#     [0.35, 1, 0],
#     [0.35, 0.35, 1]
# ]

# # Define the edges of the tetrahedron
# edges = [
#     [0, 1],
#     [0, 2],
#     [0, 3],
#     [1, 2],
#     [1, 3],
#     [2, 3]
# ]

# # Create a 3D plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Plot the vertices
# for vertex in vertices:
#     ax.scatter(vertex[0], vertex[1], vertex[2], c='r')

# # Plot the edges with colors
# for edge in edges:
#     ax.plot([vertices[edge[0]][0], vertices[edge[1]][0]],
#             [vertices[edge[0]][1], vertices[edge[1]][1]],
#             [vertices[edge[0]][2], vertices[edge[1]][2]], c='g', marker='o')

# # Set the plot limits
# ax.set_xlim([0, 1])
# ax.set_ylim([0, 1])
# ax.set_zlim([0, 1])

# # Show the plot
# plt.show()


# # Show the plot
# plt.show()



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Tetrahedron vertices
tetra_vertices = np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0], [0.5, 0.5, np.sqrt(2)/2]])


# Plot tetrahedron
c1 = (105./255.,212./255,207./255)
ax.scatter(tetra_vertices[:, 0], tetra_vertices[:, 1], tetra_vertices[:, 2], c=c1, marker='o')
c2 = (39./255.,98./255,94./255)

ax.plot_trisurf(tetra_vertices[:, 0], tetra_vertices[:, 1], tetra_vertices[:, 2], linewidth=0.5, antialiased=True, color=c2)

plt.show()
#ax.plot_trisurf(tetra_vertices[:, 0], tetra_vertices[:, 1], tetra_vertices[:, 2], linewidth=0.5, antialiased=True, color='bluegreen')
