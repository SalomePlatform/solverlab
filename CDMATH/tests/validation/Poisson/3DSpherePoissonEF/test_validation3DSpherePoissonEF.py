import FiniteElements3DPoissonSphere
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import time, json

convergence_synthesis=dict(FiniteElements3DPoissonSphere.test_desc)
def test_validation3DSphereEF():
    start = time.time()
    #### 3D sphere FE triangle mesh
    meshList=['meshSphere_1','meshSphere_2','meshSphere_3','meshSphere_4','meshSphere_5']
    meshType="Unstructured_3D_triangles"
    testColor="Green"
    nbMeshes=len(meshList)
    error_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    mesh_path='../../../ressources/3DSphere/'
    mesh_name='SphereWithTriangles'
    diag_data=[0]*nbMeshes
    resolution=100
    plt.close('all')
    i=0
    # Storing of numerical errors and mesh sizes
    for filename in meshList:
        error_tab[i], mesh_size_tab[i], min_sol_num, max_sol_num, time_tab[i] =FiniteElements3DPoissonSphere.solve(mesh_path+filename, resolution,meshType,testColor)
        assert min_sol_num>-1.1 
        assert max_sol_num<1.1
        error_tab[i]=log10(error_tab[i])
        time_tab[i]=log10(time_tab[i])
        with open('./FiniteElementsOnSpherePoisson_PlotOnSortedLines'+meshType+str(mesh_size_tab[i])+'0.csv') as f:
            lines = f.readlines()
            y = [float(line.split(",")[0]) for line in lines[1:]]
            x = [float(line.split(",")[1]) for line in lines[1:]]

        plt.plot(x, y, label= str(mesh_size_tab[i]) + ' nodes')
        mesh_size_tab[i] = 0.5*log10(mesh_size_tab[i])
        i=i+1

    end = time.time()

    # Plot over diagonal line
    plt.legend()
    plt.xlabel('Position on slice circle')
    plt.ylabel('Value on slice circle')
    plt.title('Plot over slice circle for finite elements \n for Laplace operator on 3D sphere meshes')
    plt.savefig(mesh_name+"_3DSpherePoissonFE_Slice.png")
    
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    b1=np.dot(error_tab,mesh_size_tab)   
    b2=np.sum(error_tab)
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation3DSphereEF() : Make sure you use distinct meshes and at least two meshes'
    a=( a3*b1-a2*b2)/det
    b=(-a2*b1+a1*b2)/det
    
    print( "FE on 3D sphere triangle mesh : scheme order is ", -a)
    assert abs(a+0.816)<0.1

    # Plot of convergence curves
    plt.close()
    plt.plot(mesh_size_tab, error_tab, label='log(|numerical-exact|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='least square slope : '+'%.3f' % a)
    plt.legend()
    plt.xlabel('log(sqrt(number of nodes))')
    plt.ylabel('log(error)')
    plt.title('Convergence of finite elements for \n Laplace operator on 3D sphere triangular meshes')
    plt.savefig(mesh_name+"_3DSpherePoissonFE_ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time))')
    plt.legend()
    plt.xlabel('log(sqrt(number of nodes)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite elements \n for Laplace operator on 3D sphere triangular meshes')
    plt.savefig(mesh_name+"_3DSpherePoissonFE_ComputationalTime.png")
    
    plt.close('all')

    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=[10**x for x in mesh_size_tab]
    convergence_synthesis["Space_dimension"]=3
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_cell_type"]="3DTriangles"
    convergence_synthesis["Errors"]=[10**x for x in error_tab]
    convergence_synthesis["Scheme_order"]=-a
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_Poisson_32DFV_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    test_validation3DSphereEF()
