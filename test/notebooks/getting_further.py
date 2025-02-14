#!/usr/bin/env python
# coding: utf-8

# # Complete study with Salome
# 
# ## Introduction: non-stationary 3D diffusion
# 
# In this new example, we're trying to solve an unsteady problem,
# i.e. with a temporal evolution of the solution.
# 
# The diffusion equation is again considered:
# 
# $$
# \rho c_p \frac{\partial T}{\partial t} - \lambda\Delta T = \Phi 
# $$
# 
# We begin by creating the geometry and mesh with Salome.
# 
# ## Preparing the mesh

# In[ ]:


import cdmath
import CoreFlows as cf
from math import exp
try:
    from IPython.display import Image, display
    withJupyter = True
except ImportError:
    withJupyter = False

import salome

try:
    from pvsimple import *
except:
    from paraview.simple import *


# In[ ]:


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


# In[ ]:


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


# ## Resolution with CoreFlows
# 
# Load the MED format mesh created earlier.

# In[5]:


# Mesh
mesh = cdmath.Mesh("Mesh_1.med")


# We can then define the problem with the material characteristics

# In[6]:


density = 10000
specificHeat = 300
conductivity = 5
problem = cf.DiffusionEquation(3, False, density, specificHeat, conductivity)


# specifying that a finite volume resolution is desired (argument `False`).
# We define the initial condition, detailing the definition domain,
# the list of values to be assigned (here we have only one) and finally the mesh entity type, which here is the cell
# in the finite volume framework. If we had asked for a finite element resolution, we would have had to
# give a field to the nodes as in the previous example.

# In[7]:


problem.setInitialFieldConstant(mesh, [15], cdmath.CELLS)


# As far as the boundary conditions are concerned, we choose here to give a zero flow
# on the sphere's boundaries, i.e. a Neumann condition:

# In[8]:


# Boundary conditions
problem.setNeumannBoundaryCondition("External_Face")


# We prepare a heat source centered at the sphere's origin by traversing the mesh cells
# to assign values to the field.

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


# We define some resolution parameters
# - implicit time scheme
# - matrix inversion algorithm
# - problem name
# - output frequency

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


# The resolution can then be launched:

# In[11]:


problem.initialize()
exitStatus = problem.run()
if not exitStatus:
    raise RuntimeError(f"Simulation failed")

problem.terminate()


# The results are automatically recorded according to the information
# and can be observed with **Salome**.

# In[ ]:


if withJupyter:
    diffusionEquation_diffusion_3Dpvd = PVDReader(registrationName='DiffusionEquation_diffusion_3D.pvd', FileName='./DiffusionEquation_diffusion_3D.pvd')
    animationScene1 = GetAnimationScene()
    timeKeeper1 = GetTimeKeeper()
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    renderView1 = GetActiveViewOrCreate('RenderView')
    diffusionEquation_diffusion_3DpvdDisplay = Show(diffusionEquation_diffusion_3Dpvd, renderView1, 'UnstructuredGridRepresentation')
    diffusionEquation_diffusion_3DpvdDisplay.Representation = 'Surface'
    animationScene1.GoToLast()
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
    renderView1.CameraPosition = [4.013230503985284, 1.4341159836510837, 5.1589435816460725]
    renderView1.CameraFocalPoint = [-0.00011444091796874981, -4.827976226806646e-06, -6.239227298405378e-20]
    renderView1.CameraViewUp = [-0.0835939881493395, 0.9749749758733987, -0.20599961545098572]
    renderView1.CameraParallelScale = 1.7319295439543239

    SaveScreenshot('result.png', renderView1, ImageResolution=[3056, 1115],
        TransparentBackground=1)
    display(Image(filename="result.png"))

