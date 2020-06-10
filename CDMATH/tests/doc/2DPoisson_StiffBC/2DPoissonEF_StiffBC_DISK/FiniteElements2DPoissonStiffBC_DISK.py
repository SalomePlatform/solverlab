# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Poisson 2D -\triangle u = 0 sur un disque avec conditions aux limites de Dirichlet discontinues
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Utilisation de la méthode des éléménts finis P1 avec champs u et f discrétisés aux noeuds d'un maillage triangulaire
#		        Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#               Comparaison de la solution numérique avec la solution exacte u=atan(2*x/(x*x+y*y-1))
#================================================================================================================================

import cdmath
from math import atan, pi
from numpy import sign
import time, json
import PV_routines
import VTK_routines

test_desc={}
test_desc["Initial_data"]="None"
test_desc["Boundary_conditions"]="Dirichlet"
test_desc["Global_name"]="FE simulation of the 2D Poisson equation"
test_desc["Global_comment"]="Condition limite discontinue, Maillage triangulaire"
test_desc["PDE_model"]="Poisson"
test_desc["PDE_is_stationary"]=True
test_desc["PDE_search_for_stationary_solution"]=False
test_desc["Numerical_method_name"]="P1 FE"
test_desc["Numerical_method_space_discretization"]="Finite elements"
test_desc["Numerical_method_time_discretization"]="None"
test_desc["Mesh_is_unstructured"]=True
test_desc["Geometry"]="Disk"
test_desc["Part_of_mesh_convergence_analysis"]=True

# Fonction qui remplace successivement les colonnes d'une matrices par un vecteur donné et retourne la liste des déterminants
def gradientNodal(M, values):
    matrices=[0]*(len(values)-1)
    for i in range(len(values)-1):
        matrices[i] = M.deepCopy()        
        for j in range(len(values)):
            matrices[i][j,i] = values[j]

    result = cdmath.Vector(len(values)-1)    
    for i in range(len(values)-1):
        result[i] = matrices[i].determinant()

    return result

def solve(filename,resolution,meshType, testColor):
    start = time.time()
    test_desc["Mesh_type"]=meshType
    test_desc["Test_color"]=testColor
    
    #Chargement du maillage triangulaire du disque unité D(0,1)
    #=======================================================================================
    my_mesh = cdmath.Mesh(filename+".med")
    if( my_mesh.getSpaceDimension()!=2 or my_mesh.getMeshDimension()!=2) :
        raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 2")
    if(not my_mesh.isTriangular()) :
        raise ValueError("Wrong cell types : mesh is not made of triangles")
    
    nbNodes = my_mesh.getNumberOfNodes()
    nbCells = my_mesh.getNumberOfCells()
    
    test_desc["Space_dimension"]=my_mesh.getSpaceDimension()
    test_desc["Mesh_dimension"]=my_mesh.getMeshDimension()
    test_desc["Mesh_number_of_elements"]=my_mesh.getNumberOfNodes()
    test_desc["Mesh_cell_type"]=my_mesh.getElementTypesNames()

    print("Mesh loading done")
    print("Number of nodes=", nbNodes)
    print("Number of cells=", nbCells)
    
    #Détermination des noeuds intérieurs
    #======================================================================
    my_ExactSol = cdmath.Field("Exact_field", cdmath.NODES, my_mesh, 1)
    nbInteriorNodes = 0
    nbBoundaryNodes = 0
    maxNbNeighbours = 0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
    interiorNodes=[]
    boundaryNodes=[]
    eps=1e-6
    
    #parcours des noeuds pour discrétisation du second membre et extraction 1) des noeuds intérieur 2) des noeuds frontière 3) du nb max voisins d'un noeud
    for i in range(nbNodes):
        Ni=my_mesh.getNode(i)
        x = Ni.x()
        y = Ni.y()
    
        #Robust calculation of atan(2x/(x**2+y**2-1)
        #my_ExactSol[i]=atan2(2*x*sign(x**2+y**2-1),abs(x**2+y**2-1))#mettre la solution exacte de l'edp
        if x**2+y**2-1 > eps :
            print("!!! Warning Mesh ", meshType," !!! Node is not in the unit disk.",", eps=",eps, ", x**2+y**2 - 1=",x**2+y**2 - 1)
            #raise ValueError("!!! Domain should be the unit disk.")
        if x**2+y**2-1 < -eps :
            my_ExactSol[i] = atan(2*x/(x**2+y**2-1))
        elif x>0 : #x**2+y**2-1>=0
            my_ExactSol[i] = -pi/2
        elif x<0 : #x**2+y**2-1>=0
            my_ExactSol[i] =  pi/2
        else : #x=0
            my_ExactSol[i] = 0

        if my_mesh.isBorderNode(i): # Détection des noeuds frontière
            boundaryNodes.append(i)
            nbBoundaryNodes=nbBoundaryNodes+1
        else: # Détection des noeuds intérieurs
            interiorNodes.append(i)
            nbInteriorNodes=nbInteriorNodes+1
            maxNbNeighbours= max(1+Ni.getNumberOfCells(),maxNbNeighbours) #true only in 2D, need a function Ni.getNumberOfNeighbourNodes()

    test_desc["Mesh_max_number_of_neighbours"]=maxNbNeighbours
    
    print("Right hand side discretisation done")
    print("Number of interior nodes=", nbInteriorNodes)
    print("Number of boundary nodes=", nbBoundaryNodes)
    print("Max number of neighbours=", maxNbNeighbours)
    
    # Construction de la matrice de rigidité et du vecteur second membre du système linéaire
    #=======================================================================================
    Rigidite=cdmath.SparseMatrixPetsc(nbInteriorNodes,nbInteriorNodes,maxNbNeighbours) # warning : third argument is maximum number of non zero coefficients per line of the matrix
    RHS=cdmath.Vector(nbInteriorNodes)
    
    M=cdmath.Matrix(3,3)
    # Vecteurs gradient de la fonction de forme associée à chaque noeud d'un triangle (hypothèse 2D)
    GradShapeFunc0=cdmath.Vector(2)
    GradShapeFunc1=cdmath.Vector(2)
    GradShapeFunc2=cdmath.Vector(2)
    
    #On parcourt les triangles du domaine
    for i in range(nbCells):
    
        Ci=my_mesh.getCell(i)
    
        #Contribution à la matrice de rigidité
        nodeId0=Ci.getNodeId(0)
        nodeId1=Ci.getNodeId(1)
        nodeId2=Ci.getNodeId(2)
    
        N0=my_mesh.getNode(nodeId0)
        N1=my_mesh.getNode(nodeId1)
        N2=my_mesh.getNode(nodeId2)
    
        M[0,0]=N0.x()
        M[0,1]=N0.y()
        M[0,2]=1
        M[1,0]=N1.x()
        M[1,1]=N1.y()
        M[1,2]=1
        M[2,0]=N2.x()
        M[2,1]=N2.y()
        M[2,2]=1
    
        #Values of each shape function at each node
        values0=[1,0,0]
        values1=[0,1,0]
        values2=[0,0,1]
    
        GradShapeFunc0 = gradientNodal(M,values0)/2
        GradShapeFunc1 = gradientNodal(M,values1)/2
        GradShapeFunc2 = gradientNodal(M,values2)/2
    
        #Création d'un tableau (numéro du noeud, gradient de la fonction de forme)
        GradShapeFuncs={nodeId0 : GradShapeFunc0}
        GradShapeFuncs[nodeId1]=GradShapeFunc1
        GradShapeFuncs[nodeId2]=GradShapeFunc2
    
    
        # Remplissage de  la matrice de rigidité et du second membre
        for j in [nodeId0,nodeId1,nodeId2] :
            if boundaryNodes.count(j)==0 : #seuls les noeuds intérieurs contribuent au système linéaire (matrice de rigidité et second membre)
                j_int=interiorNodes.index(j)#indice du noeud j en tant que noeud intérieur
                #Pas de contribution au second membre car pas de terme source
                boundaryContributionAdded=False#Needed in case j is a border cell
                #Contribution de la cellule triangulaire i à la ligne j_int du système linéaire
                for k in [nodeId0,nodeId1,nodeId2] : 
                    if boundaryNodes.count(k)==0 : #seuls les noeuds intérieurs contribuent à la matrice du système linéaire
                        k_int=interiorNodes.index(k)#indice du noeud k en tant que noeud intérieur
                        Rigidite.addValue(j_int,k_int,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())
                    elif boundaryContributionAdded == False: #si condition limite non nulle au bord (ou maillage non recouvrant), ajouter la contribution du bord au second membre de la cellule j
                        if boundaryNodes.count(nodeId0)!=0 :
                            u0=my_ExactSol[nodeId0]
                        else:
                            u0=0
                        if boundaryNodes.count(nodeId1)!=0 :
                            u1=my_ExactSol[nodeId1]
                        else:
                            u1=0
                        if boundaryNodes.count(nodeId2)!=0 :
                            u2=my_ExactSol[nodeId2]
                        else:
                            u2=0
                        boundaryContributionAdded=True#Contribution from the boundary to matrix line j is done in one step
                        GradGh = gradientNodal(M,[u0,u1,u2])/2
                        RHS[j_int] += -(GradGh*GradShapeFuncs[j])/Ci.getMeasure()

    print("Linear system matrix building done")
    
    # Résolution du système linéaire
    #=================================
    LS=cdmath.LinearSolver(Rigidite,RHS,100,1.E-6,"CG","LU")#Remplacer CG par CHOLESKY pour solveur direct
    LS.setComputeConditionNumber()
    SolSyst=LS.solve()
    
    print( "Preconditioner used : ", LS.getNameOfPc() )
    print( "Number of iterations used : ", LS.getNumberOfIter() )
    print("Linear system solved")
    
    test_desc["Linear_solver_algorithm"]=LS.getNameOfMethod()
    test_desc["Linear_solver_preconditioner"]=LS.getNameOfPc()
    test_desc["Linear_solver_precision"]=LS.getTolerance()
    test_desc["Linear_solver_maximum_iterations"]=LS.getNumberMaxOfIter()
    test_desc["Linear_system_max_actual_iterations_number"]=LS.getNumberOfIter()
    test_desc["Linear_system_max_actual_error"]=LS.getResidu()
    test_desc["Linear_system_max_actual_condition number"]=LS.getConditionNumber()


    # Création du champ résultat
    #===========================
    my_ResultField = cdmath.Field("ResultField", cdmath.NODES, my_mesh, 1)
    for j in range(nbInteriorNodes):
        my_ResultField[interiorNodes[j]]=SolSyst[j];#remplissage des valeurs pour les noeuds intérieurs
    for j in range(nbBoundaryNodes):
        my_ResultField[boundaryNodes[j]]=my_ExactSol[boundaryNodes[j]];#remplissage des valeurs pour les noeuds frontière (condition limite)
    #sauvegarde sur le disque dur du résultat dans un fichier paraview
    my_ResultField.writeVTK("FiniteElements2DPoissonStiffBC_DISK_"+meshType+str(nbNodes))
    my_ExactSol.writeVTK("ExactSol2DPoissonStiffBC_DISK_"+meshType+str(nbNodes))
    
    print("Numerical solution of 2D Poisson equation on a disk using finite elements done")
    
    end = time.time()

    #Calcul de l'erreur commise par rapport à la solution exacte
    #===========================================================
    l2_norm_sol_exacte=my_ExactSol.normL2()[0]
    l2_error = (my_ExactSol - my_ResultField).normL2()[0]
    
    print("L2 relative error = norm( exact solution - numerical solution )/norm( exact solution ) = ", l2_error/l2_norm_sol_exacte)
    print("Maximum numerical solution = ", my_ResultField.max(), " Minimum numerical solution = ", my_ResultField.min())
    print("Maximum exact solution = ", my_ExactSol.max(), " Minimum exact solution = ", my_ExactSol.min())
    
    #Postprocessing : 
    #================
    # Extraction of the diagonal data
    diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,-1,0],[0,1,0], resolution)
    # save 2D picture
    PV_routines.Save_PV_data_to_picture_file("FiniteElements2DPoissonStiffBC_DISK_"+meshType+str(nbNodes)+'_0.vtu',"ResultField",'NODES',"FiniteElements2DPoissonStiffBC_DISK_"+meshType+str(nbNodes))
    PV_routines.Save_PV_data_to_picture_file("ExactSol2DPoissonStiffBC_DISK_"+meshType+str(nbNodes)+'_0.vtu',"Exact_field",'NODES',"ExactSol2DPoissonStiffBC_DISK_"+meshType+str(nbNodes))
    
    test_desc["Computational_time_taken_by_run"]=end-start
    test_desc["Absolute_error"]=l2_error
    test_desc["Relative_error"]=l2_error/l2_norm_sol_exacte

    with open('test_PoissonStiffBC'+str(my_mesh.getMeshDimension())+'D_EF_'+"DISK_"+meshType+str(nbCells)+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)

    return l2_error/l2_norm_sol_exacte, nbNodes, diag_data, my_ResultField.min(), my_ResultField.max(), end - start

if __name__ == """__main__""":
    solve("diskWithTriangles",100,"Unstructured_triangles","Green")
