"""
Remark that in Strain2RigidCosseratMapping the default elongation is 0 
but for the Frames2StrainCosseratMapping the default elongation is fixed to 1. 

"""
import numpy as np
import matplotlib.pyplot as plt

data_frames2strain = np.loadtxt("monitor_frames2strain_x.txt")
data_strain2rigid = np.loadtxt("monitor_strain2rigid_x.txt")

time = data_strain2rigid[:, 0]
strains_f2s = data_frames2strain[:, 1:] # for Frames2Strain
strains_s2r = data_strain2rigid[:, 1:] # for Strain2Rigid


# On applique la correction sur toutes les sections
nb_sections = 8


nb_iterations = len(time)
print(f"Nombre d'itérations : {nb_iterations}")
# 1. Erreur globale
erreur_global = np.linalg.norm(strains_f2s - strains_s2r, axis=1)

plt.figure(figsize=(10, 6))
plt.plot(time, erreur_global, label="Erreur globale (norme L2)", color="red")
# plt.yscale("log")
plt.xlabel("Time")
plt.ylabel("Erreur")
plt.title("Erreur entre Frames2Strain et Strain2Rigid")
plt.legend()
plt.grid()
plt.savefig(f"erreur_global{nb_sections}_setprecision.png", dpi=300)
plt.show()


# 2. Erreur par section (à implémenter)

f2s_3d = strains_f2s.reshape(nb_iterations, nb_sections, 6)
s2r_3d = strains_s2r.reshape(nb_iterations, nb_sections, 6)

erreur_par_section = np.linalg.norm(f2s_3d - s2r_3d, axis=2)

# for i in range(nb_sections):
#     plt.figure(figsize=(10, 6))
#     plt.plot(np.arange(nb_iterations), erreur_par_section[:, i], label=f"Erreur section {i}", color="blue")
#     plt.xlabel("Time")
#     plt.ylabel("Erreur")
#     plt.title(f"Erreur entre Frames2Strain et Strain2Rigid - Section {i}")
#     plt.legend()
#     plt.grid()
#     plt.show()

max_erreur_par_section = np.max(erreur_par_section, axis=0)
plt.figure(figsize=(10, 6))
plt.bar(range(nb_sections), max_erreur_par_section, color="blue")
plt.xlabel("Section")
plt.ylabel("Max Erreur")
plt.title("Max Erreur entre Frames2Strain et Strain2Rigid par section")
plt.xticks(range(nb_sections))
plt.grid()
plt.savefig(f"max_erreur_par_section{nb_sections}_setprecision.png", dpi=300)
plt.show()