//============================================================================
/**
 * \file TransportEquation.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date 24 March 2015
 * \brief Fluid enthalpy transport equation
 * */
//============================================================================

/*! \class TransportEquation TransportEquation.hxx "TransportEquation.hxx"
 *  \brief Scalar advection equation for a fluid enthalpy
 *  \details see \ref TransportEqPage for more details
 */
#ifndef TransportEquation_HXX_
#define TransportEquation_HXX_

#include "ProblemCoreFlows.hxx"

using namespace std;

class TransportEquation: public ProblemCoreFlows
{

public :
	/** \fn TransportEquation
			 * \brief Constructor for the enthalpy transport in a fluid
			 * \param [in] phaseType : \ref Liquid or \ref Gas
			 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
			 * \param [in] vector<double> : fluid velocity (assumed constant)
			 *  */
	TransportEquation(phaseType fluid, pressureEstimate pEstimate,vector<double> vitesseTransport);

	//Gestion du calcul
	virtual void initialize();
	virtual void terminate();//vide la mémoire et enregistre le résultat final
	bool initTimeStep(double dt);
	double computeTimeStep(bool & stop);//propose un pas de temps pour le calcul. Celà nécessite de discrétiser les opérateur (convection, diffusion, sources) et pour chacun d'employer la condition cfl. En cas de problème durant ce calcul (exemple t=tmax), renvoie stop=true
	void abortTimeStep();//efface les inconnues calculées par solveTimeStep() et reinitialise dt à 0
	bool iterateTimeStep(bool &ok);
	virtual void save();
	virtual void validateTimeStep();

	/** \fn setIntletBoundaryCondition
			 * \brief adds a new boundary condition of type Inlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [out] void
			 *  */
	void setInletBoundaryCondition(string groupName,double enthalpy){
		_limitField[groupName]=LimitField(Inlet,-1,vector<double>(_Ndim,0),vector<double>(_Ndim,0),vector<double>(_Ndim,0),-1,enthalpy,-1,-1);
	};

	/*Physical parameters*/
	Field& getFluidTemperatureField(){
		return _TT;
	}

	void setLiqSatEnthalpy(double hsatl){
		_hsatl=hsatl;
	};
	void setVapSatEnthalpy(double hsatv){
		_hsatv=hsatv;
	};
	void setLiqSatDensity(double rhosatl){
		_rhosatl=rhosatl;
	};
	void setVapSatDensity(double rhosatv){
		_rhosatv=rhosatv;
	};
	void setTransportVelocity(Vector v){
		_vitesseTransport=v;
	};
protected :
	double computeTransportMatrix();
	double computeRHS();
	void updatePrimitives();
	double temperature(double h){
		return _Tref+(h-_href)/_cpref;
	};
	double voidFraction(double h){
		double titre=(h-_href)/(_hsatv-_hsatl);
		if (titre<0)
			return 0;
		else if (titre>1)
			return 1;
		else return titre*_rhosatl/(titre*_rhosatl+(1-titre)*_rhosatv);
	};
	double density(double alpha){
		return alpha*_rhosatv+(1-alpha)*_rhosatl;
	};

	Field   _TT, _Alpha, _Rho;//Fields of temperature and coupled temperature
	double _rhosatv, _rhosatl;
	double _Tref, _href, _cpref;
	Vector _vitesseTransport, _normale;
	bool _transportMatrixSet;
	Vec _Hn, _deltaH, _Hk, _Hkm1, _b0;
	double _dt_transport, _dt_src;
};

#endif /* TransportEquation_HXX_ */
