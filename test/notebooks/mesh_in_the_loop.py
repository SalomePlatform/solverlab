#!/usr/bin/env python
# coding: utf-8

# # Mesh in the loop
# 
# This exercise demonstrates how to use the Solverlab module to refine and solve a problem in several steps.

# In[ ]:


from pathlib import Path
from math import sin, pi, tanh, cos
import sys
try:
    from IPython.display import Image, display
    withJupyter = True
except ImportError:
    withJupyter = False

original_stdout = sys.stdout

import salome
# salome.salome_init()
# from salome.shaper import model
# import SHAPERSTUDY
# import SMESH, SALOMEDS
# from salome.smesh import smeshBuilder

import CoreFlows as cf
import cdmath

try:
    import pvsimple
    from pvsimple import *
except:
    import paraview.simple as pvsimple
    from paraview.simple import *

from salome.shaper import model
import SHAPERSTUDY
import SMESH, SALOMEDS
from salome.smesh import smeshBuilder

sys.stdout = original_stdout


# ## Problem definition
# 
# We are looking to solve a 2D Poisson problem on a domain:
# 
# $$
# \forall x \in \Omega, \ \Delta u(x) = -f(x)\\
# \forall x \in \partial \Omega, \ u(x) = u_d(x)
# $$
# 
# where $u$ is the unknown, $f$ the source term and $u_d$ the Dirichlet boundary conditions.
# 
# To make it easier to write the solution to the problem, the terms of the equation are calculated in the following function. The terms of this particular problem are obtained from the following manufactured solution:
# $$
# u_{m}(x) = \mathrm{tanh}\left(-100 \left(y - \frac{1}{2} - \frac{1}{4} \mathrm{sin}\left(2 \pi x\right) \right) \right) + \mathrm{tanh} \left(100 \left(y - x \right) \right)
# $$
# 
# ![Manufacture solution](manufactured_solution.png)

# In[87]:


def equation_terms(x, y):
    A = tanh(-100 * (y - 0.5 - 0.25 * sin(2 * pi * x)))
    B = tanh(100 * (y - x))

    dAx = 50 * pi * cos(2 * pi * x) * (1 - A**2)
    dAxx = 100 * pi * (pi * sin(2 * pi * x) * (A**2 - 1) - cos(2 * pi * x) * dAx * A)
    dAy = 100 * (A**2 - 1)
    dAyy = 200 * dAy * A

    dBx = 100 * (B**2 - 1)
    dBxx = 200 * dBx * B
    dBy = 100 * (1 - B**2)
    dByy = -200 * dBy * B

    fm = dAxx + dAyy + dBxx + dByy
    um = A + B

    return (-um, -fm)


# ## Initial mesh
# 
# The first step is to create a very coarse mesh. We use the **Shaper** module to create the geometry, then the Smesh module to automatically create the mesh.
# 
# The geometry is a simple square of size 1.

# In[88]:


model.begin()
partSet = model.moduleDocument()
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))
SketchLine_1 = Sketch_1.addLine(1, 0, 0, 0)
SketchProjection_1 = Sketch_1.addProjection(
    model.selection("VERTEX", "PartSet/Origin"), False
)
SketchPoint_1 = SketchProjection_1.createdFeature()
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchPoint_1.result(), True)
SketchLine_2 = Sketch_1.addLine(0, 0, 0, 1.000000000001039)
SketchLine_3 = Sketch_1.addLine(0, 1.000000000001039, 1, 1)
SketchLine_4 = Sketch_1.addLine(1, 1, 1, 0)

Sketch_1.setCoincident(SketchLine_4.endPoint(), SketchLine_1.startPoint(), True)
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchLine_2.startPoint(), True)
Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchLine_3.startPoint(), True)
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_4.startPoint(), True)
Sketch_1.setPerpendicular(SketchLine_1.result(), SketchLine_2.result(), True)
Sketch_1.setPerpendicular(SketchLine_2.result(), SketchLine_3.result(), True)
Sketch_1.setPerpendicular(SketchLine_3.result(), SketchLine_4.result(), True)
Sketch_1.setHorizontal(SketchLine_1.result(), True)
Sketch_1.setLength(SketchLine_4.result(), 1, True)
Sketch_1.setLength(SketchLine_1.result(), 1, True)

model.do()
Face_1 = model.addFace(
    Part_1_doc,
    [
        model.selection(
            "FACE",
            "Sketch_1/Face-SketchLine_4r-SketchLine_3r-SketchLine_2r-SketchLine_1r",
        )
    ],
)
model.end()

model.publishToShaperStudy()
(Face_1_1,) = SHAPERSTUDY.shape(model.featureStringId(Face_1))


# We now create the mesh from the geometry. We create an initial, relatively coarse mesh, with a deliberately small number of segments per side of the square (`Nsegments = 50`).
# 
# The mesh can be exported in MED format for visualization. This mesh is the starting point for a resolution of satisfactory quality in several steps.

# In[89]:


Nsegments = 50
smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Face_1_1, "Mesh_1")
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Simple_Parameters_1 = NETGEN_1D_2D.Parameters(smeshBuilder.SIMPLE)
NETGEN_2D_Simple_Parameters_1.SetNumberOfSegments(Nsegments)
NETGEN_2D_Simple_Parameters_1.LengthFromEdges()
NETGEN_2D_Simple_Parameters_1.SetAllowQuadrangles(0)
Mesh_1.Compute()
Mesh_1.CheckCompute()
smesh.SetName(Mesh_1, "Mesh_1")

initial_mesh_file = Path("initial_mesh.med")
Mesh_1.ExportMED(str(initial_mesh_file), 0, 41, 1, Mesh_1, 1, [], "", -1, 1)


# ## Remesh in the loop
# 
# ### Initial mesh resolution
# 
# In this section, we solve the problem on the relatively coarse mesh we've just built.
# To do this, we use the **CoreFlows** module (from **Salome/Solverlab**). A small function is prepared so that we can easily restart the resolution with other meshes later.

# In[ ]:


def solve_poisson_problem_with_coreflows(mesh_file: Path):
    
    mesh = cdmath.Mesh(str(mesh_file))
    boundaryIds = mesh.getBoundaryNodeIds()
    dirichletValues = {}
    for i in boundaryIds:
        Ni = mesh.getNode(i)
        x = Ni.x()
        y = Ni.y()
        (u, f) = equation_terms(x, y)
        dirichletValues[i] = u

    problem = cf.StationaryDiffusionEquation(2, True, 1)
    problem.setMesh(mesh)
    problem.setDirichletValues(cf.MapIntDouble(dirichletValues))

    sourceField = cdmath.Field("sourceField", cdmath.NODES, mesh, 1)
    for i in range(mesh.getNumberOfNodes()):
        Ni = mesh.getNode(i)
        x = Ni.x()
        y = Ni.y()

        (u, f) = equation_terms(x, y)
        sourceField[i] = -f
    
    problem.setFluidTemperature(0)
    problem.setHeatPowerField(sourceField)
    problem.setFileName(mesh_file.stem)
    problem.setLinearSolver(cf.CG, cf.CHOLESKY)
    problem.initialize()
    problem.solveStationaryProblem()
    
    simSolution = problem.getOutputTemperatureField()
    simSolution.writeMED(f"{mesh_file.stem}_simulationSolution")

solve_poisson_problem_with_coreflows(initial_mesh_file)


# ### Show solution
# 
# The results can be observed using the Paravis module. Similarly, we write a function to be able to observe other results later.

# In[ ]:


def export_solution_png(solution_file: Path):
    pvsimple._DisableFirstRenderCameraReset()
    
    solution = MEDReader(
        registrationName='solution',
        FileNames=[str(solution_file.resolve())]
    )
    
    renderView1 = GetActiveViewOrCreate('RenderView')
    display = Show(solution, renderView1, 'UnstructuredGridRepresentation')
    display.Representation = 'Surface'
    renderView1.ResetCamera(False)
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5, 0.5, 10000.0]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
    materialLibrary1 = GetMaterialLibrary()
    renderView1.Update()

    ColorBy(display, ('POINTS', 'Temperature'))
    display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    temperatureLUT = GetColorTransferFunction('Temperature')
    temperaturePWF = GetOpacityTransferFunction('Temperature')
    temperatureTF2D = GetTransferFunction2D('Temperature')

    # change representation type
    display.SetRepresentationType('Surface With Edges')
    renderView1.OrientationAxesVisibility = 0
    display.SetScalarBarVisibility(renderView1, True)
    temperatureLUTColorBar = GetScalarBar(temperatureLUT, renderView1)
    temperatureLUT.RescaleTransferFunction(-2.0, 2.0)
    
    temperatureLUTColorBar.TitleFontSize = 8
    temperatureLUTColorBar.LabelFontSize = 8
    temperatureLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
    temperatureLUTColorBar.LabelColor = [0.0, 0.0, 0.0]

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5, 0.5, 10000.0]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
    renderView1.CameraParallelScale = 0.7071067811865476

    # save screenshot
    image_file = Path(solution_file.stem + ".png")
    SaveScreenshot(str(image_file), renderView1, ImageResolution=[1368, 1020], TransparentBackground=1)
    

first_solve_result_file = Path(initial_mesh_file.stem +"_simulationSolution.med")


# In[ ]:


if withJupyter:
    export_solution_png(first_solve_result_file)
    display(Image(filename=first_solve_result_file.stem+".png"))


# The resulting solution doesn't really resemble the manufactured solution. In this case, we can try to refine the mesh to observe a more precise solution. You can recreate a finer mesh as at the beginning of this exercise, but perhaps try using the first solution obtained to refine the mesh only where relevant.
# 
# ### Mesh size field
# 
# Before remeshing, we need to define a new mesh size field. To do this, we'll use the resulting solution. As you can see, the particularity of the problem stems from the very pronounced gradient between the different uniform zones of the manufactured solution. We therefore need to calculate a mesh size that is large in the uniform zones and small in the high-gradient zones.
# 
# In the following, we use the following formula to calculate the mesh size field.
# 
# $$
# h_{adapted}(x) = \frac{h_{max}}{ \left\vert \left\vert \nabla u(x) \right\vert \right\vert + 1 }
# $$
# 
# where $h_{max}$ is the maximum mesh size and $a$ is a parameter such that $h_{min} = h_{max}/a$ and $h_{min}$ is the minimum desired mesh cell size.

# In[93]:


def compute_mesh_size_field(solution_file: Path, mesh_size_file: Path, h_max, a=20):
    from paraview.simple import MEDReader,Gradient, PointDatatoCellData, Calculator, SaveData
    
    field1 = MEDReader(registrationName="field1", FileNames=[str(solution_file.resolve())])
    gradient1 = Gradient(
        registrationName="gradient1", Input=field1, ScalarArray="Temperature"
    )
    cellGradient1 = PointDatatoCellData(
        registrationName="cellGradient1", Input=gradient1
    )
    calculator1 = Calculator(registrationName="Calculator1", Input=cellGradient1)
    calculator1.AttributeType = "Cell Data"
    calculator1.ResultArrayName = "h_adapted"
    calculator1.Function = f"{h_max}/min(sqrt(Gradient_X^2+Gradient_Y^2)+1, {a})"

    SaveData(str(mesh_size_file.resolve()), proxy=calculator1)


# Here, we choose the parameter $h_{max} = 0.3$ because we want to keep a relatively coarse mesh where the solution is uniform and we arbitrarily choose the parameter $a=30$ to have a minimum mesh size 30 times smaller than the maximum size.

# In[94]:


second_mesh_file = Path("second_mesh.med")
second_mesh_size_field = Path(second_mesh_file.stem + "_size_field.med")
compute_mesh_size_field(first_solve_result_file, second_mesh_size_field, 0.3, a=30)


# ### Remeshing from the initial mesh
# 
# Now, to remesh, simply call up the MGAdapt function provided in the **Smesh** module.

# In[95]:


def remesh_mgadapt(mesh_size_file: Path, new_mesh_file: Path):
    smesh = smeshBuilder.New()
    objet_adapt = smesh.Adaptation("MG_Adapt")
    objet_adapt.setMEDFileIn(str(mesh_size_file.resolve()))
    objet_adapt.setMEDFileOut(str(new_mesh_file.resolve()))

    hypo = smesh.CreateAdaptationHypothesis()
    hypo.setSizeMapType("Background")
    objet_adapt.setMEDFileBackground(str(mesh_size_file.resolve()))
    hypo.setSizeMapFieldName("h_adapted")
    hypo.setTimeStepRank(0, 0)

    objet_adapt.AddHypothesis(hypo)
    objet_adapt.Compute(False)

    # resave med in 2D space
    ([Mesh_1], status) = smesh.CreateMeshesFromMED(str(new_mesh_file))
    smesh.SetName(Mesh_1, 'Mesh_1')
    Mesh_1.ExportMED(str(new_mesh_file), 0, 41, 1, Mesh_1, 1, [], '',-1, 1 )
    

remesh_mgadapt(second_mesh_size_field, second_mesh_file)


# ### Resolution on refined mesh
# 
# The problem can be solved again with the new mesh.

# In[ ]:


solve_poisson_problem_with_coreflows(second_mesh_file)
second_solve_result_file = Path(second_mesh_file.stem +"_simulationSolution.med")


# In[ ]:


if withJupyter:
    export_solution_png(second_solve_result_file)
    display(Image(filename=second_solve_result_file.stem+".png"))


# The new mesh gives much better results. We can have fun trying a new mesh based on the second solution. Since the second solution is better, the gradient estimate for calculating the mesh size field should also be more accurate. We can therefore use more ambitious remeshing parameters.

# In[ ]:


third_mesh_file = Path("third_mesh.med")
third_mesh_size_field = Path(third_mesh_file.stem + "_size_field.med")
compute_mesh_size_field(second_solve_result_file, third_mesh_size_field, 0.2, a=60)

remesh_mgadapt(third_mesh_size_field, third_mesh_file)

solve_poisson_problem_with_coreflows(third_mesh_file)
third_solve_result_file = Path(third_mesh_file.stem +"_simulationSolution.med")


# In[ ]:


if withJupyter:
    export_solution_png(third_solve_result_file)
    display(Image(filename=third_solve_result_file.stem+".png"))

