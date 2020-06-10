/*
 * fieldtests.cxx
 *
 *  Created on: 24 janv. 2012
 *      Authors: CDMAT
 */

#include "FieldTests.hxx"
#include "Vector.hxx"
#include <string>

#include <MEDCouplingFieldDouble.hxx>

using namespace std;
using namespace MEDCoupling;

//----------------------------------------------------------------------
void
FieldTests::testClassField( void )
//----------------------------------------------------------------------
{
	Mesh M(0.0,1.0,10,0.,1.,5);

	Field conc1("CONCENTRATION",CELLS,M,2,1.2) ;
	CPPUNIT_ASSERT_EQUAL( 1.2, conc1.getTime() );
    for (int j=0;j<conc1.getNumberOfComponents();j++)
    	for (int i=0;i<conc1.getNumberOfElements();i++)
    			conc1(i,j)=i+j;
    string fileNameVTK="champ";
    conc1.writeVTK(fileNameVTK);

    string fileNameMED="champ";
    conc1.writeMED(fileNameMED);
    conc1.setTime(2.3,1);
    conc1.writeMED(fileNameMED,false);
    for (int j=0;j<conc1.getNumberOfComponents();j++)
    	for (int i=0;i<conc1.getNumberOfElements();i++)
    		CPPUNIT_ASSERT_EQUAL( double(i+j), conc1(i,j) );
	CPPUNIT_ASSERT_EQUAL( 2, conc1.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc1.getNumberOfElements() );
	CPPUNIT_ASSERT_EQUAL( 2.3, conc1.getTime() );

	Field conc1n("CONCENTRATION",NODES,M,2,1.2) ;
	CPPUNIT_ASSERT_EQUAL( 1.2, conc1n.getTime() );
    for (int j=0;j<conc1n.getNumberOfComponents();j++)
    	for (int i=0;i<conc1n.getNumberOfElements();i++)
    		conc1n(i,j)=i*1.0;
    string fileNameVTKn="champn";
    conc1n.writeVTK(fileNameVTKn);

    string fileNameMEDn="champn";
    conc1n.writeMED(fileNameMEDn);
    conc1n.setTime(2.3,1);
    conc1n.writeMED(fileNameMEDn,false);

    for (int j=0;j<conc1n.getNumberOfComponents();j++)
    	for (int i=0;i<conc1n.getNumberOfElements();i++)
    		CPPUNIT_ASSERT_EQUAL( 1.0*i, conc1n(i,j) );
	CPPUNIT_ASSERT_EQUAL( 2, conc1n.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 66, conc1n.getNumberOfElements() );
	CPPUNIT_ASSERT_EQUAL( 2.3, conc1n.getTime() );

	Field conc6("CONCENTRATION",CELLS,M,2);
    for (int i=0;i<conc6.getNumberOfComponents();i++)
    	for (int j=0;j<conc6.getNumberOfElements();j++)
    		conc6(j,i)=i*1.0+2.*j;

    for (int i=0;i<conc6.getNumberOfComponents();i++)
        for (int j=0;j<conc6.getNumberOfElements();j++)
        	CPPUNIT_ASSERT_EQUAL( 1.0*i+2.*j, conc6.getValues()[i+j*conc6.getNumberOfComponents()] );

    CPPUNIT_ASSERT_EQUAL( 2, conc6.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc6.getNumberOfElements() );

	Field conc6n("CONCENTRATION",NODES,M,2);
    for (int i=0;i<conc6n.getNumberOfComponents();i++)
    	for (int j=0;j<conc6n.getNumberOfElements();j++)
    		conc6n(j,i)=i*1.0+2.*j;

    for (int i=0;i<conc6n.getNumberOfComponents();i++)
        for (int j=0;j<conc6n.getNumberOfElements();j++)
        	CPPUNIT_ASSERT_EQUAL( 1.0*i+2.*j, conc6n.getValues()[i+j*conc6n.getNumberOfComponents()] );

    CPPUNIT_ASSERT_EQUAL( 2, conc6n.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 66, conc6n.getNumberOfElements() );

	Field conc3(conc1) ;
    for (int j=0;j<conc6.getNumberOfComponents();j++)
    	for (int i=0;i<conc3.getNumberOfElements();i++)
    		conc3(i,j)=-(i+j);

	Vector v1=conc3.getValuesOnComponent(1);
	Vector v2=conc3.getValuesOnAllComponents(4);

	for (int i=0;i<conc3.getNumberOfElements();i++)
		CPPUNIT_ASSERT_EQUAL( double(-(i+1)), v1(i) );

	for (int j=0;j<conc3.getNumberOfComponents();j++)
		CPPUNIT_ASSERT_EQUAL( double(-(4+j)), v2(j) );

	double x=conc3(2,0);
	CPPUNIT_ASSERT_EQUAL( x, -2.0 );

	for (int i=0;i<conc3.getNumberOfElements();i++)
		CPPUNIT_ASSERT_EQUAL( double(-i), conc3(i) );
	CPPUNIT_ASSERT_EQUAL( 2, conc3.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc3.getNumberOfElements() );

	conc6=conc3+conc1;
    for (int i=0;i<conc6.getNumberOfElements();i++)
    	CPPUNIT_ASSERT_EQUAL( 0., conc6[i] );
	CPPUNIT_ASSERT_EQUAL( 2, conc6.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc6.getNumberOfElements() );

	conc6=conc3-conc1;
    for (int i=0;i<conc6.getNumberOfElements();i++)
    	CPPUNIT_ASSERT_EQUAL( -2.*(i), conc6(i) );
	CPPUNIT_ASSERT_EQUAL( 2, conc6.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc6.getNumberOfElements() );

	conc6=conc1;
	conc6+=conc1;
    for (int j=0;j<conc6.getNumberOfComponents();j++)
    	for (int i=0;i<conc6.getNumberOfElements();i++)
    		CPPUNIT_ASSERT_EQUAL( 2.0*(i+j), conc6(i,j) );
	CPPUNIT_ASSERT_EQUAL( 2, conc6.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc6.getNumberOfElements() );

	conc6=conc1;
	conc6*=2.0;
    for (int j=0;j<conc6.getNumberOfComponents();j++)
    	for (int i=0;i<conc6.getNumberOfElements();i++)
    		CPPUNIT_ASSERT_EQUAL( 2.0*(i+j), conc6(i,j) );
	CPPUNIT_ASSERT_EQUAL( 2, conc6.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc6.getNumberOfElements() );

	Field conc7("CONCENTRATION",CELLS,M,2) ;
	MCAuto<MEDCouplingFieldDouble> f1=conc1.getField();
	conc7.setFieldByMEDCouplingFieldDouble(f1);
    conc7.setName("CONC");
    for (int i=0;i<conc7.getNumberOfElements();i++)
    {
    	CPPUNIT_ASSERT_EQUAL( 1.0*i, conc7(i) );
    	CPPUNIT_ASSERT_EQUAL( 1.0*i, conc7[i] );
    }
	CPPUNIT_ASSERT_EQUAL( 2, conc7.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc7.getNumberOfElements() );
	CPPUNIT_ASSERT( conc7.getName().compare("CONC")==0 );

	Field conc77("CONCENTRATION",CELLS,M,2) ;
	conc77.setInfoOnComponent(0,"compo1");
	conc77.setInfoOnComponent(1,"compo2");
	CPPUNIT_ASSERT(conc77.getInfoOnComponent(0).compare("compo1")==0 );
	CPPUNIT_ASSERT(conc77.getInfoOnComponent(1).compare("compo2")==0 );

	MCAuto<MEDCouplingFieldDouble> f2=conc1.getField();
	conc77.setFieldByDataArrayDouble(f2->getArray());
    conc77.setName("CONC");
    for (int i=0;i<conc77.getNumberOfElements();i++)
    {
    	CPPUNIT_ASSERT_EQUAL( 1.0*i, conc77(i) );
    	CPPUNIT_ASSERT_EQUAL( 1.0*i, conc77[i] );
    }
	CPPUNIT_ASSERT_EQUAL( 2, conc77.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc77.getNumberOfElements() );
	CPPUNIT_ASSERT( conc77.getName().compare("CONC")==0 );

	Field conc8("CONCENTRATION",CELLS,M) ;
    for (int i=0;i<conc8.getNumberOfElements();i++)
    	conc8[i]=i*1.0;
    for (int i=0;i<conc8.getNumberOfElements();i++)
    {
    	CPPUNIT_ASSERT_EQUAL( 1.0*i, conc8(i) );
    	CPPUNIT_ASSERT_EQUAL( 1.0*i, conc8[i] );
    }
	CPPUNIT_ASSERT_EQUAL( 1, conc8.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 50, conc8.getNumberOfElements() );

	Field conc8n("CONCENTRATION",NODES,M) ;
    for (int i=0;i<conc8n.getNumberOfElements();i++)
    	conc8n[i]=i*1.0;
    for (int i=0;i<conc8n.getNumberOfElements();i++)
    {
    	CPPUNIT_ASSERT_EQUAL( 1.0*i, conc8n(i) );
    	CPPUNIT_ASSERT_EQUAL( 1.0*i, conc8n[i] );
    }
	CPPUNIT_ASSERT_EQUAL( 1, conc8n.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 66, conc8n.getNumberOfElements() );

	Field conc9=conc8 ;
	conc9/=2.0;
    for (int i=0;i<conc9.getNumberOfElements();i++)
    	CPPUNIT_ASSERT_EQUAL( 1.0*i/2., conc9(i) );

    Field conc10 ;
	conc10=conc8;
	conc10-=2.0;
    for (int i=0;i<conc10.getNumberOfElements();i++)
    	CPPUNIT_ASSERT_EQUAL( 1.0*i-2.0, conc10(i) );

    Field conc11=conc8 ;
	conc11+=2.0;
    for (int i=0;i<conc11.getNumberOfElements();i++)
    	CPPUNIT_ASSERT_EQUAL( 1.0*i+2.0, conc11(i) );

    Field conc12=conc8 ;
	conc12+=conc8;
    for (int i=0;i<conc12.getNumberOfElements();i++)
    	CPPUNIT_ASSERT_EQUAL( 2.0*i, conc12(i) );

    Field conc13=conc8 ;
	conc13-=conc8;
    for (int i=0;i<conc13.getNumberOfElements();i++)
    	CPPUNIT_ASSERT_EQUAL( 0.0, conc13(i) );

    Field conc14=2.*conc1 ;
    Field conc15=conc1*2. ;
    Field conc16=conc1/3. ;

    for (int i=0;i<conc14.getNumberOfElements();i++)
    {
    	CPPUNIT_ASSERT_EQUAL( conc1(i)*2., conc14(i) );
    	CPPUNIT_ASSERT_EQUAL( conc1(i)*2., conc15(i) );
    	CPPUNIT_ASSERT_EQUAL( conc1(i)/3., conc16(i) );
    }

	Mesh MF(0.0,1.0,3,0.,1.,3);
	Field concF1("CONCENTRATION",FACES,MF) ;
    for (int j=0;j<concF1.getNumberOfComponents();j++)
    	for (int i=0;i<concF1.getNumberOfElements();i++)
    		concF1(i,j)=i+j;

    for (int j=0;j<concF1.getNumberOfComponents();j++)
    	for (int i=0;i<concF1.getNumberOfElements();i++)
    		CPPUNIT_ASSERT_EQUAL( double(i+j), concF1(i,j) );
	CPPUNIT_ASSERT_EQUAL( 1, concF1.getNumberOfComponents() );
	CPPUNIT_ASSERT_EQUAL( 0.0, concF1.getTime() );
	CPPUNIT_ASSERT_EQUAL( 24, concF1.getNumberOfElements() );
}
