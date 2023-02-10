#!/usr/bin/env python3
# -*-coding:utf-8 -*

import sys
import solverlab
from math import exp

#===============================================================================================================================
# Name        : Simulation of a 3D heat equation 
# Description : Test solving the diffusion of the temperature T in a solid from a non constant initial datum
#               \rho cp dT/dt-\lambda\Delta T=\Phi 
#               Dirichlet boundary conditions
#               Finite elements or finite volumes
#               Initial temperature is stronger at the center of the ball
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2022
#================================================================================================================================


def DiffusionEquation_3D_BALL_Dirichlet(FECalculation, fileName):

	""" Description : Test solving the diffusion of the temperature T in a solid from a non constant initial datum. 
		Equation : Thermal diffusion equation  \rho cp dT/dt-\lambda\Delta T=\Phi.
		        Heat capacity cp, density \rho, and conductivity \lambda of the solid MUST be defined.
                Initial temperature is stronger at the center of the ball
	"""
	#Space dimension of the problem
	spaceDim=3
	
	# Prepare for the mesh
	print("Loading mesh " );
	ball_mesh = solverlab.Mesh("../resources/Ball_1.med")
	if( ball_mesh.getSpaceDimension()!=3 or ball_mesh.getMeshDimension()!=3) :
	    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 3")
	if(not ball_mesh.isTetrahedral()) :
		raise ValueError("Wrong cell types : mesh is not made of tetrahedra")
	
	nbNodes = ball_mesh.getNumberOfNodes()
	nbCells = ball_mesh.getNumberOfCells()

	print( "Loaded a 3D ball mesh with ", nbNodes, " nodes and ",nbCells, " cells")

	# Definition of field support parameter
	if( FECalculation):
		supportOfField=solverlab.NODES
	else:
		supportOfField=solverlab.CELLS	

	boundaryTemperature = 15
	#Set the rinitial temperature
	initialField = solverlab.Field("Initial temperature", supportOfField, ball_mesh, 1)
	if( FECalculation):
		for i in range(ball_mesh.getNumberOfNodes()):
			Ni= ball_mesh.getNode(i)
			x = Ni.x()
			y = Ni.y()
			z = Ni.z()
	
			initialField[i]=boundaryTemperature+50*exp(-(x*x+y*y+z*z))#mettre la fonction definie au second membre de l'edp
	else:
		for i in range(ball_mesh.getNumberOfCells()):
			Ci= ball_mesh.getCell(i)
			x = Ci.x()
			y = Ci.y()
			z = Ci.z()
	
			initialField[i]=boundaryTemperature+50*exp(-(x*x+y*y+z*z))#mettre la fonction definie au second membre de l'edp
	initialField.writeVTK("InitialTemperature")
		
    # Mandatory physical values
	solid_specific_heat=300# specific heat capacity, default value 300
	solid_density=10000# density, default value 10000
	solid_conductivity=5# conductivity, default value 5

	myProblem = solverlab.DiffusionEquation(spaceDim,FECalculation,solid_density,solid_specific_heat,solid_conductivity);
	myProblem.setInitialField(initialField);

	# set the Dirichlet boundary condition
	myProblem.setDirichletBoundaryCondition("Sphere",boundaryTemperature)
	
    # Set the mesh and initial data
	#initial_data_inputfile="../resources/Ball_1";
	#initial_data_fieldName="Solid temperature";
	#print("Loading unstructured mesh and initial data", " in file ", initial_data_inputfile )
	#initial_data_time_iteration=0# default value is 0
	#initial_data_time_sub_iteration=0# default value is 0
	#initial_data_time_meshLevel=0# default value is 0
	#myProblem.setInitialField(initial_data_inputfile, initial_data_fieldName, initial_data_time_iteration, initial_data_time_sub_iteration, initial_data_time_meshLevel, supportOfField)


    # set the numerical method
	myProblem.setTimeScheme( solverlab.Implicit)
	max_nb_its_lin_solver = 50
	myProblem.setLinearSolver(solverlab.GMRES, solverlab.LU, max_nb_its_lin_solver );

    # computation parameters
	MaxNbOfTimeStep = 3 ;# default value is 10
	freqSave = 1;# default value is 1
	cfl = 100;# default value is 1
	maxTime = 100000000;# default value is 10
	precision = 1e-6;# default value is 1e-6
	result_directory="."# default value = current directory

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setResultDirectory(result_directory)
	#myProblem.setSaveFileFormat(solverlab.MED)#default value is solverlab.VTK

    # evolution
	myProblem.initialize();
	print("Running python test "+ fileName );

	ok = myProblem.run();
	if (ok):
		print( "Python simulation " + fileName + " is successful !" );
		pass
	else:
		print( "Python simulation " + fileName + "  failed ! " );
		pass

	print( "------------ End of simulation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        FECalculation = bool(int(sys.argv[1]))
        # name of result file
        if( FECalculation):
        	fileName = "3D_BALL_FE_Dirichlet";# default value is ""
        else:
        	fileName = "3D_BALL_FV_Dirichlet";# default value is ""
        DiffusionEquation_3D_BALL_Dirichlet(FECalculation, fileName)
    else :
        raise ValueError("DiffusionEquation_3D_BALL_Dirichlet : missing one argument")
