import FiniteVolumes3DPoisson_CUBE
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json

convergence_synthesis=dict(FiniteVolumes3DPoisson_CUBE.test_desc)

def test_validation3DVF_checkerboard():
    start = time.time()
    #### 3D FV cartesian non conforming mesh
    meshList=['cubeWithCheckerboardCubes_2x2x2','cubeWithCheckerboardCubes_4x4x4','cubeWithCheckerboardCubes_8x8x8','cubeWithCheckerboardCubes_16x16x16']
    meshType="Non_conforming_checkerboard"
    testColor="Orange, BC violated. PB with mesh ?"
    nbMeshes=len(meshList)
    error_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    mesh_path='./'
    mesh_name='CubeWithCheckerboardCubes'
    diag_data=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(3),resolution+1)
    plt.close('all')
    i=0
    
    print( "!!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!!!")
    print( "Mesh was badly generated (look for instance  at boundary groups)")
    
    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
        error_tab[i], mesh_size_tab[i], diag_data[i], min_sol_num, max_sol_num, time_tab[i] =FiniteVolumes3DPoisson_CUBE.solve_file(mesh_path+filename,resolution,meshType,testColor)
        assert min_sol_num>-1.1 
        assert max_sol_num<8.1
        plt.plot(curv_abs, diag_data[i], label= str(mesh_size_tab[i]) + ' cells')
        error_tab[i]=log10(error_tab[i])
        time_tab[i]=log10(time_tab[i])
        i=i+1
        
    end = time.time()

    # Plot over diagonal line
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Value on diagonal line')
    plt.title('Plot over diagonal line for finite volumes \n for Laplace operator on 3D checkerboard meshes')
    plt.savefig(mesh_name+"_3DPoissonVF_PlotOverDiagonalLine.png")

    for i in range(nbMeshes):
        mesh_size_tab[i] = 1./3*log10(mesh_size_tab[i])

    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    b1=np.dot(error_tab,mesh_size_tab)   
    b2=np.sum(error_tab)
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation3DVF_checkerboard() : Make sure you use distinct meshes and at least two meshes'
    a=( a3*b1-a2*b2)/det
    b=(-a2*b1+a1*b2)/det
    
    print( "FV on 3D checkerboard mesh : scheme order is ", -a)
    assert abs(a-0.26)<0.01
    
    # Plot of convergence curve
    plt.close()
    plt.plot(mesh_size_tab, error_tab, label='log(|numerical-exact|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='least square slope : '+'%.3f' % a)
    plt.legend()
    plt.xlabel('1/3 log(number of cells)')
    plt.ylabel('log(error)')
    plt.title('Convergence of finite volumes for \n Laplace operator on 3D checkerboard meshes')
    plt.savefig(mesh_name+"_3DPoissonVF_ConvergenceCurve.png")

    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('1/3 log(number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes \n for Laplace operator on 3D checkerboard meshes')
    plt.savefig(mesh_name+"_3DPoissonVF_ComputationalTime.png")
    
    plt.close('all')

    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=[10**x for x in mesh_size_tab]
    convergence_synthesis["Space_dimension"]=3
    convergence_synthesis["Mesh_dimension"]=3
    convergence_synthesis["Mesh_cell_type"]="Cubes"
    convergence_synthesis["Errors"]=[10**x for x in error_tab]
    convergence_synthesis["Scheme_order"]=-a
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_Poisson_3DFV_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    test_validation3DVF_checkerboard()
