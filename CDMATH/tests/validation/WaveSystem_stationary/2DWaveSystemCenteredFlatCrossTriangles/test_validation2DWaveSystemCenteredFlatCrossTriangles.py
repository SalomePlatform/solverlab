import WaveSystemCentered
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import sys
import time, json

    
def test_validation2DWaveSystemCenteredFlatCrossTriangles(scaling):
    start = time.time()
    #### 2D flat cross triangle meshes
    meshList=['squareWithFlatCrossTriangles_00','squareWithFlatCrossTriangles_0','squareWithFlatCrossTriangles_1']#,'squareWithFlatCrossTriangles_2','squareWithFlatCrossTriangles_3'
    meshType="Unstructured_flat_cross_triangles"
    testColor="Green"
    nbMeshes=len(meshList)
    mesh_size_tab=[0]*nbMeshes
    mesh_path='../../../ressources/2DFlatCrossTriangles/'
    mesh_name='squareWithFlatCrossTriangles'
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)

    error_p_tab=[0]*nbMeshes
    error_u_tab=[0]*nbMeshes
    diag_data_press=[0]*nbMeshes
    diag_data_vel=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    t_final=[0]*nbMeshes
    ndt_final=[0]*nbMeshes
    max_vel=[0]*nbMeshes
    cond_number=[0]*nbMeshes

    plt.close('all')
    i=0
    cfl=0.5
    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
        error_p_tab[i], error_u_tab[i], mesh_size_tab[i], t_final[i], ndt_final[i], max_vel[i], diag_data_press[i], diag_data_vel[i], time_tab[i], cond_number[i] = WaveSystemCentered.solve_file(mesh_path+filename, mesh_name, resolution,scaling,meshType,testColor,cfl,"Periodic")
        assert max_vel[i]>0.94 and max_vel[i]<1.1
        error_p_tab[i]=log10(error_p_tab[i])
        error_u_tab[i]=log10(error_u_tab[i])
        time_tab[i]=log10(time_tab[i])
        i=i+1
    
    end = time.time()

    # Plot over diagonal line
    for i in range(nbMeshes):
        if(scaling==0):
            plt.plot(curv_abs, diag_data_press[i], label= str(mesh_size_tab[i]) + ' cells - no scaling')
        else:
            plt.plot(curv_abs, diag_data_press[i], label= str(mesh_size_tab[i]) + ' cells - with scaling')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Pressure on diagonal line')
    plt.title('Plot over diagonal line for stationary wave system \n with centered scheme on 2D flat cross triangular meshes')
    plt.savefig(mesh_name+'_Pressure_2DWaveSystemCentered_'+"scaling"+str(scaling)+"_PlotOverDiagonalLine.png")
    plt.close()

    plt.clf()
    for i in range(nbMeshes):
        if(scaling==0):
            plt.plot(curv_abs, diag_data_vel[i],   label= str(mesh_size_tab[i]) + ' cells - no scaling')
        else:
            plt.plot(curv_abs, diag_data_vel[i],   label= str(mesh_size_tab[i]) + ' cells - with scaling')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Velocity on diagonal line')
    plt.title('Plot over diagonal line for the stationary wave system \n with centered scheme on 2D flat cross triangle meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_PlotOverDiagonalLine.png")    
    plt.close()

    # Plot of number of time steps
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, ndt_final, label='Number of time step to reach stationary regime - no scaling')
    else:
        plt.plot(mesh_size_tab, ndt_final, label='Number of time step to reach stationary regime - with scaling')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Max time steps for stationary regime')
    plt.title('Number of times steps required for the stationary Wave System \n with centered scheme on 2D flat cross triangle meshes')
    plt.savefig(mesh_name+"_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_TimeSteps.png")
    
    # Plot of number of stationary time
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, t_final, label='Time where stationary regime is reached - no scaling')
    else:
        plt.plot(mesh_size_tab, t_final, label='Time where stationary regime is reached - with scaling')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Max time for stationary regime')
    plt.title('Simulated time for the stationary Wave System \n with centered scheme on 2D flat cross triangle meshes')
    plt.savefig(mesh_name+"_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_FinalTime.png")
    
    # Plot of number of maximal velocity norm
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, max_vel, label='Maximum velocity norm - no scaling')
    else:
        plt.plot(mesh_size_tab, max_vel, label='Maximum velocity norm - with scaling')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Max velocity norm')
    plt.title('Maximum velocity norm for the stationary Wave System \n with centered scheme on 2D flat cross triangle meshes')
    plt.savefig(mesh_name+"_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_MaxVelNorm.png")
    
    # Plot of condition number 
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, cond_number, label='Condition number - no scaling')
    else:
        plt.plot(mesh_size_tab, cond_number, label='Condition number - with scaling')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Condition number')
    plt.title('Condition number for the stationary Wave System \n with centered scheme on 2D square meshes')
    plt.savefig(mesh_name+"_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_condition_number.png")
    
    for i in range(nbMeshes):
        mesh_size_tab[i]=0.5*log10(mesh_size_tab[i])
        
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DWaveSystemFVCentered() : Make sure you use distinct meshes and at least two meshes'

    b1p=np.dot(error_p_tab,mesh_size_tab)   
    b2p=np.sum(error_p_tab)
    ap=( a3*b1p-a2*b2p)/det
    bp=(-a2*b1p+a1*b2p)/det
    
    if(scaling==0):
        print("FV Centered on 2D flat cross triangle meshes : scheme order for pressure without scaling is ", -ap)
    else:
        print("FV Centered on 2D flat cross triangle meshes : scheme order for pressure with scaling is ", -ap)

    b1u=np.dot(error_u_tab,mesh_size_tab)   
    b2u=np.sum(error_u_tab)
    au=( a3*b1u-a2*b2u)/det
    bu=(-a2*b1u+a1*b2u)/det
    
    if(scaling==0):
        print("FV Centered on 2D flat cross triangle meshes : scheme order for velocity without scaling is ", -au)
    else:
        print("FV Centered on 2D flat cross triangle meshes : scheme order for velocity with scaling is ", -au)
    
    # Plot of convergence curves
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, error_p_tab, label='|error on stationary pressure| - no scaling')
    else:
        plt.plot(mesh_size_tab, error_p_tab, label='|error on stationary pressure| - with scaling')
    plt.legend()
    plt.xlabel('1/2 log(Number of cells)')
    plt.ylabel('log(|error p|)')
    plt.title('Convergence of finite volumes for the stationary Wave System \n with centered scheme on 2D flat cross triangle meshes')
    plt.savefig(mesh_name+"_Pressure_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_ConvergenceCurve.png")
    
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, error_u_tab, label='log(|error on stationary velocity|) - no scaling')
    else:
        plt.plot(mesh_size_tab, error_u_tab, label='log(|error on stationary velocity|) - with scaling')
    plt.legend()
    plt.xlabel('1/2 log(Number of cells)')
    plt.ylabel('log(|error u|)')
    plt.title('Convergence of finite volumes for the stationary Wave System \n with centered scheme on 2D flat cross triangle meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, time_tab, label='log(cpu time) - no scaling')
    else:
        plt.plot(mesh_size_tab, time_tab, label='log(cpu time) - with scaling')
    plt.legend()
    plt.xlabel('1/2 log(Number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for the stationary Wave System \n with centered scheme on 2D flat cross triangle meshes')
    plt.savefig(mesh_name+"2DWaveSystemCentered_"+"scaling"+str(scaling)+"_ComputationalTime.png")

    plt.close('all')

    convergence_synthesis={}

    convergence_synthesis["PDE_model"]="Wave system"
    convergence_synthesis["PDE_is_stationary"]=False
    convergence_synthesis["PDE_search_for_stationary_solution"]=True
    if(scaling==0):
        convergence_synthesis["Numerical_method_name"]="Centered scheme no scaling"
    else:
        convergence_synthesis["Numerical_method_name"]="Centered scheme scaling"			
    convergence_synthesis["Numerical_method_space_discretization"]="Finite volumes"
    convergence_synthesis["Numerical_method_time_discretization"]="Implicit"
    convergence_synthesis["Initial_data"]="Constant pressure, divergence free velocity"
    convergence_synthesis["Boundary_conditions"]="Periodic"
    convergence_synthesis["Numerical_parameter_cfl"]=cfl
    convergence_synthesis["Space_dimension"]=2
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=mesh_size_tab
    convergence_synthesis["Mesh_cell_type"]="Triangles"
    convergence_synthesis["Numerical_error_velocity"]=error_u_tab
    convergence_synthesis["Numerical_error_pressure"]=error_p_tab
    convergence_synthesis["Max_vel_norm"]=max_vel
    convergence_synthesis["Final_time"]=t_final  
    convergence_synthesis["Final_time_step"]=ndt_final  
    convergence_synthesis["Scheme_order"]=-au
    convergence_synthesis["Scheme_order_vel"]=-au
    convergence_synthesis["Scaling_preconditioner"]=scaling
    convergence_synthesis["Condition_numbers"]=cond_number
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_WaveSystem_2DFV_Centered_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        scaling = int(sys.argv[1])
        test_validation2DWaveSystemCenteredFlatCrossTriangles(scaling)
    else :
        test_validation2DWaveSystemCenteredFlatCrossTriangles(2)

