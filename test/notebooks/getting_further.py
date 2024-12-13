#!/usr/bin/env python
# coding: utf-8

# # Étude complète avec Salome
# 
# ## Introduction: diffusion instationnaire en 3D
# 
# Dans ce nouvel exemple, on cherche à résoudre un problème instationnaire,
# c'est-à-dire avec une évolution temporelle de la solution.
# 
# On considère à nouveau l'équation de diffusion :
# 
# $$
# \rho c_p \frac{\partial T}{\partial t} - \lambda\Delta T = \Phi 
# $$
# 
# On commence par créer la géométrie et le maillage avec Salome.
# 
# ## Préparation du maillage

# In[1]:


import cdmath
import CoreFlows as cf
from math import exp

import salome
#salome.salome_init()

#import pvsimple
#pvsimple.ShowParaviewView()
from pvsimple import *
#pvsimple._DisableFirstRenderCameraReset()


# In[2]:


# Create the geometry
from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

Sphere_1 = model.addSphere(Part_1_doc, model.selection("VERTEX", "PartSet/Origin"), 1)

model.end()
model.publishToShaperStudy()

import SHAPERSTUDY
Sphere_1_1, = SHAPERSTUDY.shape(model.featureStringId(Sphere_1))


# In[3]:


# Create the mesh from the geometry
import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

Mesh_1 = smesh.Mesh(Sphere_1_1,'Mesh_1')
GMSH_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH_3D)
Gmsh_3D_Parameters = GMSH_3D.Parameters()
Gmsh_3D_Parameters.Set3DAlgo( 0 )
Gmsh_3D_Parameters.SetMinSize( 0.05 )
Gmsh_3D_Parameters.SetIs2d( 0 )
status = Mesh_1.RemoveHypothesis(Gmsh_3D_Parameters)
Gmsh_3D_Parameters_1 = GMSH_3D.Parameters()
Gmsh_3D_Parameters_1.Set3DAlgo( 0 )
Gmsh_3D_Parameters_1.SetMinSize( 0 )
Gmsh_3D_Parameters_1.SetSizeFactor( 0.05 )
Gmsh_3D_Parameters_1.SetIs2d( 0 )
GMSH_2D = Mesh_1.Triangle(algo=smeshBuilder.GMSH_2D)
Gmsh_Parameters = GMSH_2D.Parameters()
Gmsh_Parameters.Set2DAlgo( 0 )
Gmsh_Parameters.SetMinSize( 0 )
Gmsh_Parameters.SetMaxSize( 0.1)
Gmsh_Parameters.SetIs2d( 1 )
isDone = Mesh_1.Compute()
Mesh_1.CheckCompute()


# In[4]:


Group_1 = Mesh_1.CreateEmptyGroup( SMESH.FACE, 'External_Face' )
nbAdd = Group_1.AddFrom( Mesh_1.GetMesh() )
[ Group_1 ] = Mesh_1.GetGroups()

smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Gmsh_3D_Parameters, 'Gmsh 3D Parameters')
smesh.SetName(GMSH_3D.GetAlgorithm(), 'GMSH_3D')
smesh.SetName(Gmsh_3D_Parameters_1, 'Gmsh 3D Parameters')
smesh.SetName(GMSH_2D.GetAlgorithm(), 'GMSH_2D')

try:
  Mesh_1.ExportMED( r'Mesh_1.med', 0, 41, 1, Mesh_1, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')


# ## Résolution avec CoreFlows
# 
# On charge le maillage au format MED créé juste avant.

# In[5]:


# Mesh
mesh = cdmath.Mesh("Mesh_1.med")


# On peut alors définir le problème avec les caractéristiques du matériau

# In[6]:


density = 10000
specificHeat = 300
conductivity = 5
problem = cf.DiffusionEquation(3, False, density, specificHeat, conductivity)


# en précisant qu'on souhaite ici une résolution par les volume finis (argument `False`).
# On définit la condition initiale en détaillant le domaine de définition,
# la liste des valeurs à assigner (ici on en a qu'une) et enfin le support, qui est ici la cellule
# dans le cadre des volumes finis. Si on avait demandé une résolution aux éléments finis, on aurait du
# donner une champ aux noeuds comme dans l'exemple précédent.

# In[7]:


problem.setInitialFieldConstant(mesh, [15], cdmath.CELLS)


# En ce qui concerne les conditions aux limites, on choisit ici de donner une flux nul
# sur les frontières de la sphère, soit une condition de Neumann :

# In[8]:


# Boundary conditions
problem.setNeumannBoundaryCondition("External_Face")


# On prépare une source de chaleur centrée à l'origine de la sphère en pacourant les
# cellules du maillage pour attribuer les valeurs au champ.

# In[9]:


# Right-hand side
heatPowerField = cdmath.Field("Heat power", cdmath.CELLS, mesh, 1)
for i in range(mesh.getNumberOfCells()):
    Ci = mesh.getCell(i)
    x = Ci.x()
    y = Ci.y()
    z = Ci.z()

    heatPowerField[i] = 1000 * exp(-10 * (x * x + y * y + z * z))

problem.setHeatPowerField(heatPowerField)
problem.initialize()
exitStatus = problem.run()
if not exitStatus:
    raise RuntimeError(f"Simulation failed")

problem.terminate()


# On définit quelques paramètres de résolution
# - schéma temporel implicite
# - l'algorithme d'inversion de matrice
# - le nom du problème
# - la fréquence de sortie

# In[10]:


# Resolution
problem.setTimeScheme(cf.Implicit)
problem.setLinearSolver(cf.GMRES, cf.LU, 50)

problem.setCFL(100.0)
problem.setPrecision(1e-6)
problem.setMaxNbOfTimeStep(50)
problem.setTimeMax(100000000000)
problem.setFreqSave(100)
problem.setFileName("diffusion_3D")
problem.setResultDirectory(".")


# 
# Et on peut alors lancer la résolution :
# 

# In[11]:


problem.initialize()
exitStatus = problem.run()
if not exitStatus:
    raise RuntimeError(f"Simulation failed")

problem.terminate()


# Les résultats sont automatiquement enregistrés selon les informations données
# au problème et peuvent être observés avec **Salome**.

# In[12]:




diffusionEquation_diffusion_3Dpvd = PVDReader(registrationName='DiffusionEquation_diffusion_3D.pvd', FileName='./DiffusionEquation_diffusion_3D.pvd')


# In[13]:


animationScene1 = GetAnimationScene()
timeKeeper1 = GetTimeKeeper()
animationScene1.UpdateAnimationUsingDataTimeSteps()
renderView1 = GetActiveViewOrCreate('RenderView')
diffusionEquation_diffusion_3DpvdDisplay = Show(diffusionEquation_diffusion_3Dpvd, renderView1, 'UnstructuredGridRepresentation')
diffusionEquation_diffusion_3DpvdDisplay.Representation = 'Surface'
animationScene1.GoToLast()


# In[14]:


renderView1.ResetCamera(False)
materialLibrary1 = GetMaterialLibrary()
renderView1.Update()
clip1 = Clip(registrationName='Clip1', Input=diffusionEquation_diffusion_3Dpvd)
clip1.ClipType.Origin = [0.0, 0.0, 0.0]
clip1.ClipType.Normal = [0.0, 0.0, 1.0]
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')
clip1Display.Representation = 'Surface'
Hide(diffusionEquation_diffusion_3Dpvd, renderView1)
renderView1.Update()
ColorBy(clip1Display, ('CELLS', 'SOLVERLAB results'))
clip1Display.RescaleTransferFunctionToDataRange(True, False)
clip1Display.SetScalarBarVisibility(renderView1, True)
sOLVERLABresultsLUT = GetColorTransferFunction('SOLVERLABresults')
sOLVERLABresultsPWF = GetOpacityTransferFunction('SOLVERLABresults')
sOLVERLABresultsTF2D = GetTransferFunction2D('SOLVERLABresults')
clip1Display.RescaleTransferFunctionToDataRange(False, True)
SetActiveSource(diffusionEquation_diffusion_3Dpvd)
HideInteractiveWidgets(proxy=clip1.ClipType)
SetActiveSource(clip1)
ShowInteractiveWidgets(proxy=clip1.ClipType)
layout1 = GetLayout()
layout1.SetSize(3056, 1115)


# In[15]:


renderView1.CameraPosition = [4.013230503985284, 1.4341159836510837, 5.1589435816460725]
renderView1.CameraFocalPoint = [-0.00011444091796874981, -4.827976226806646e-06, -6.239227298405378e-20]
renderView1.CameraViewUp = [-0.0835939881493395, 0.9749749758733987, -0.20599961545098572]
renderView1.CameraParallelScale = 1.7319295439543239

SaveScreenshot('result.png', renderView1, ImageResolution=[3056, 1115],
    TransparentBackground=1)


# [Résultat diffusion 3D](result.png)

# In[ ]:




