#!/usr/bin/env python
# coding: utf-8

# 
# # Documentation SolverLab pour Salome
# 
# ## Introduction
# 
# **SOLVERLAB** se décompose en deux paquets Python:
# 
# - `cdmath` rassemble les outils élémentaires pour construire et résoudre des systèmes matriciels
# - `CoreFlows` assemble les outils pour résoudre efficacement des problèmes d'ingénierie classiques
# 
# **SOLVERLAB** se veut rassembler les librairies `MEDCOUPLING` et `PETSc` afin de constituer un outil puissant, pratique et facile à prendre en main pour résoudre des équations aux dérivées partielles.
# 
# Le projet est initialement hébergé sur [GitHub](https://github.com/ndjinga/SOLVERLAB) et une version stable est incorporée dans **Salome**.
# 
# ## Installation
# 
# **SOLVERLAB** est installé avec **Salome**. Un ensemble de ressource est disponible dans l'installation : `path-to-salome/INSTALL/SOLVERLAB/share`, notamment la documentation du module `CoreFlows` sous forme de pdf dans `SOLVERLAB/share/doc` ou sous forme `doxygen` rassemblant l'ensemble des options disponibles pour chaques modèles est disponible dans `SOLVERLAB/share/doc/coreflows-dev-doc/html/index.html`.
# 
# ## Pour commencer
# 
# Une fois **Salome** et **SOLVERLAB** installés, il suffit de charger l'environnement de **SOLVERLAB** à partir du fichier présent dans l'installation :
# 
# ```bash
# source path-to-salome/INSTALL/SOLVERLAB/env_SOLVERLAB.sh
# ```
# 
# * * *
# 
# ## Exemple: équation de diffusion en 1D
# 
# L'exemple déroulé ici repris de l'exemple d'équation de la chaleur sur [GitHub](https://github.com/ndjinga/SOLVERLAB/blob/master/CDMATH/tests/doc/1DHeatEquation/HeatEquation1D_RegularGrid.ipynb)
# 
# ### Développement en volumes finis
# 
# Dans cet exemple, on considère l'équation de diffusion suivante :
# 
# $$
# \frac{\partial u}{\partial t} = d \frac{\partial^2 u}{\partial x^2}
# $$
# 
# sur le domaine 1D, $\Omega = [0,1]$ et où $d$ est le coefficient de diffusion. On considère des conditions aux limites périodiques de chaque coté et une condition initiale de la forme $u(x,t=0) = u_0(x)$
# 
# On cherche à résoudre numériquement cette équation à l'aide de la méthode des volumes finis. Pour ce faire,  
# on commence par décomposer le domaine en N intervalles $[x_i, x_{i+1}]_{i=1..N}$ de taille $\Delta x = 1 / N$.  
# De même, le domaine en temps est discrétisé et l'intervalle de temps est noté $\Delta t_n = t_{n+1} - t_n$.  
# On cherche donc les valeurs moyennes :
# 
# $$
# u_i^n = \frac{1}{\Delta x} \int_{x_i}^{x_{i+1}}  u(x,t_n) \mathrm{d}x
# $$
# 
# sur chaque intervalle d'espace.
# 
# Le développement en volumes finis donne finalement l'équation:
# 
# $$
# \frac{u_i^{n+1} - u_i^n}{\Delta t_n} = d \frac{u_{i+1}^{n,n+1} - 2u_i^{n,n+1} + 2 u_{i-1}^{n,n+1}}{\Delta x^2}
# $$
# 
# avec $u_i^{n,n+1} = \frac{1}{\Delta t} \int_{t_n}^{t_{n+1}} u_i(t) \mathrm{d}t$. Cette quantité peut être exprimée simplement à l'aide de la formule du rectangle  
# de deux manières:
# 
# - $u_i^{n,n+1} = \frac{1}{\Delta t_n} \int_{t_n}^{t_{n+1}} u_i(t) \mathrm{d}t = u_i^n$ menant à une résolution explicite avec un critère de stabilité $\Delta t_n \lt \Delta x^2 / 2d$
# - $u_i^{n,n+1} = \frac{1}{\Delta t_n} \int_{t_n}^{t_{n+1}} u_i(t) \mathrm{d}t = u_i^{n+1}$ menant à une résolution implicite nécessitant l'inversion d'une matrice (inconditionnellement stable).
# 
# Pour la résolution explicite, on obtient :
# 
# $$
# u_i^{n+1} = u_i^{n} + \frac{d \Delta t_n}{\Delta x^2} \left( u_{i+1}^n - 2 u_{i}^n + u_{i-1}^n\right)
# $$
# 
# et pour la résolution implicite, on obtient :
# 
# $$
# \left(1 + 2 \frac{d \Delta t_n}{\Delta x^2} \right) u_i^{n+1} -  \frac{d \Delta t_n}{\Delta x^2} \left[ u_{i+1}^{n+1} + u_{i-1}^{n+1}\right] = u_i^n
# $$
# 
# ### Résolution numérique
# 
# Dans cet exemple, on va utiliser seulement le paquet `cdmath`. En premier lieu, on définit les paramètres du problème ainsi que la condition initiale `U_initial` et les coordonnées des cellules :
# 
# 

# In[1]:


from math import sin, pi
import matplotlib.pyplot as plt

import cdmath

# Parameters
N = 100  # number of elements
d = 0.1  # diffusion coefficient
dt = 0.00001  # time interval
time_steps = 100

dx = 1.0 / N
x = [0.5 * dx + i * dx for i in range(N)]
U_initial = [
    0.5 * (1 + sin(4 * pi * xi - pi * 0.5)) * int(xi < 0.5) * int(0 < xi)
    + int(0.6 < xi) * int(xi < 0.85)
    for xi in x
]


# 
# Pour la résolution explicite, la matrice s'écrit telle que $\mathbf{U}^{n+1} = \mathbf{M} \mathbf{U}^n$ :
# 
# 

# In[2]:


# Explicit resolution
print(f"CFL condition: {dt:2.3e} < {dx * dx / 2 / d:2.3e}")
mat_explicit = cdmath.SparseMatrixPetsc(N, N, 3)
for i in range(N):
    mat_explicit.setValue(i, (i + 1) % N, d * dt / dx**2)
    mat_explicit.setValue(i, i, 1 - (2 * d * dt / dx**2))
    mat_explicit.setValue(i, (i + 1) % N, d * dt / dx**2)


# 
# Il convient éventuellement de s'assurer que la condition de stabilité est bien vérifiée. Ensuite, il suffit de faire une boucle sur les pas de temps
# 
# 

# In[3]:


U_explicit = cdmath.Vector(N)
for i in range(N):
    U_explicit[i] = U_initial[i]

for nn in range(time_steps):
    U_explicit = mat_explicit * U_explicit


# 
# Pour lé résolution implicite, on écrit la matrice telle que $\mathbf{A} \mathbf{U}^{n+1} = \mathbf{U}^n$.  
# Le probème demande alors l'inversion de la matrice $\mathbf{A}$ à chaque pas de temps. On utilise pour cela l'outil de résolution de système linéaire de `cdmath`. La méthode de résolution se base sur les sous-espaces de Krylov.
# 
# 

# In[4]:


# Implicit resolution
mat_implicit = cdmath.SparseMatrixPetsc(N, N, 3)

for i in range(N):
    mat_implicit.setValue(i, (i + 1) % N, -d * dt / dx**2)
    mat_implicit.setValue(i, i, 1.0 + (2 * d * dt / dx**2))
    mat_implicit.setValue(i, (i + 1) % N, -d * dt / dx**2)

U_implicit = cdmath.Vector(N)
for i in range(N):
    U_implicit[i] = U_initial[i]

LS = cdmath.LinearSolver(
    mat_implicit, U_implicit, 50, 1e-5, "GMRES", "ILU"
)  # 50: max number of iterations, 1e-5: precision
LS.setComputeConditionNumber()

for nn in range(time_steps):
    LS.setSndMember(U_implicit)
    U_implicit = LS.solve()
    if not LS.getStatus():
        raise RuntimeError("Implicit resolution did not converge")


# 
# Enfin, on peut tracer les profiles avec `matplotlib` :
# 
# 

# In[5]:


# Post processing
fig, ax = plt.subplots()
ax.plot(x, U_initial, label="Initial")
ax.plot(x, U_explicit.getArray(), label="Explicit")
ax.plot(x, U_implicit.getArray(), linestyle="--", label="Implicit")

ax.legend(loc="upper right")
ax.set_title(f"Results at time $t={time_steps*dt:2.3e}$ s")
ax.set_xlabel("$x$")
ax.set_ylabel("$u_i$")
ax.set_xlim(0, 1)

fig.savefig("diffusion_1D.svg")


# 
# 
# ![Exemple de résultat](diffusion_1D.svg)
# 
# * * *
# 
# ## Exemple: diffusion stationnaire en 2D
# 
# Cet exemple cherche à résoudre l'équation de diffusion dans un contexte 2D:
# 
# $$
# -\lambda\Delta T=\Phi(T) + \lambda_{sf} (T_{fluid}-T)
# $$
# 
# avec
# 
# - $T$ l'inconnue, définie pour le solide
# - $\lambda$ la conductivité du solide
# - $\Phi(T)$ le terme source
# - $T_{fluid}$ la température d'un fluide et $\lambda_{sf}$ calibrant les échanges de température entre le fluide et le solide.
# 
# Le paquet `CoreFlows` permet de paramétrer très simplement la résolution de cette équation. Pour cela, on commence par définir le domaine de résolution en construisant un maillage à l'aide de l'outil de construction de `cdmath` :
# 
# 

# In[6]:


import CoreFlows as cf
import cdmath
from math import sin, pi

xinf, xsup = -0.5, 0.5
yinf, ysup = -0.5, 0.5
nx, ny = 5, 5
mesh = cdmath.Mesh(xinf, xsup, nx, yinf, ysup, ny, 0)  # Regular triangular mesh


# 
# Ce maillage est un rectangle (pavé en 3D) composé de cellules à 3 noeuds (4 en 3D). Il faut noter qu'il est aussi possible d'importer des maillages au format `MED` par exemple.
# 
# Avec ce maillage, on peut créer le *problème* qui va servir à résoudre l'équation de diffusion :
# 
# 

# In[7]:


problem = cf.StationaryDiffusionEquation(2, True, 1.75)
problem.setMesh(mesh)


# 
# La création de ce problème prend trois arguments :
# 
# - la dimmension de l'espace, ici `2`
# - un booléan indiquant si la résolution doit être faite avec la méthode des éléments finis (`True`) ou la méthode des volumes finis (`False`)
# - la valeur de la conductivité du matériau, ici `1.75`
# 
# En ce qui concerne les conditions aux limites, on doit d'abord rassembler les noeuds des bords du domaines en des entitées nommées. Dans le cas présent, on cherche à récupérer les noeuds sur un plan en donnant la valeur de la coordonnée (ex: `xsup`), le numéro de la direction (`0` pour $x$), une valeur de précision et enfin, le nom de l'entitée :
# 
# 

# In[8]:


eps = 1e-6
mesh.setGroupAtPlan(xsup, 0, eps, "Right")
mesh.setGroupAtPlan(xinf, 0, eps, "Left")
mesh.setGroupAtPlan(ysup, 1, eps, "Front")
mesh.setGroupAtPlan(yinf, 1, eps, "Back")


# 
# Puis on peut assigner les conditions aux limites :
# 
# 

# In[9]:


problem.setDirichletBoundaryCondition("Right", 0)
problem.setDirichletBoundaryCondition("Left", 0)
problem.setDirichletBoundaryCondition("Front", 0)
problem.setDirichletBoundaryCondition("Back", 0)


# 
# Ici il s'agit de conditions aux limites de Dirichlet mais il est aussi possible d'appliquer d'autres types de conditions aux limites (Neumann, ...). Il est conseillé de se référer à la documentation dans ce cas.  
# 
# Le terme source est un champ appliqué communément au second membre de l'équation (*RHS* ou *right hand side*). Ici, on voudrait appliquer un terme source dont l'intensité dépend de sa position dans le domaine. Pour cela, on commence par créer un champ en lui donnant un nom, un type (noeuds, cellules), un maillage de référence et son nombre de composante :
# 
# 

# In[10]:


source = cdmath.Field("RHS_field", cdmath.NODES, mesh, 1)


# 
# Les valeurs du champ peuvent être assignées en parcourant la liste des noeuds de son maillage:
# 
# 

# In[11]:


for i in range(mesh.getNumberOfNodes()):
    Ni = mesh.getNode(i)
    x = Ni.x()
    y = Ni.y()

    source[i] = 2 * pi * pi * sin(pi * (x - xinf)) * sin(pi * (y - yinf))


# 
# Enfin, il ne reste plus qu'a donner ce champ au problème :
# 
# 

# In[12]:


problem.setHeatPowerField(source)


# 
# En ce qui concerne le fluide, on peut définir sa température en donnant une valeur uniforme ou un champ :
# 

# In[13]:


problem.setFluidTemperature(0.5)


# 
# Avant de passer à la résolution, on peut donner un nom au problème. Ce nom sert à l'enregistrement des fichiers de résultat.
# 
# 

# In[14]:


fileName = "StationaryDiffusion_2DEF_StructuredTriangles"
problem.setFileName(fileName)


# 
# Dans le cas présent, le problème est résolution en une seule étape. On peut préciser l'algorithme de résolution ainsi que le préconditionneur à donner à `PETSc`
# 
# 

# In[15]:


problem.setLinearSolver(cf.GMRES, cf.ILU)


# 
# On peut enfin passer à la résolution :
# 
# 

# In[16]:


problem.initialize()
exitStatus = problem.solveStationaryProblem()
if not exitStatus:
    raise RuntimeError(f"{fileName} simulation failed")


# 
# Les champs de résultat (ici l'inconnue `T` appelée `Temprature` dans ce problème) sont automatiquement sauvegardé et peuvent être visionnés avec **Salome**. Il est aussi possible de récupérer le champ solution en python afin d'en extraire des propriétés intéressantes comme la valeur maximale par exemple :
# 
# 

# In[17]:


T_field = problem.getOutputTemperatureField()
print(T_field.max())


# 
# Pour terminer, le problème peut libérer la mémoire utilisée pour stocker les différents objets `PETSc` :
# 
# 

# In[18]:


problem.terminate()


# In[ ]:




