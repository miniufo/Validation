//
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class AAMDefinition{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.csm");
		CsmDescriptor csd=(CsmDescriptor)df.getDataDescriptor();
		CylindricalSpatialModel csm=new CylindricalSpatialModel(csd);
		SphericalSpatialModel ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		
		Range range=new Range("",csd);
		
		Variable u=new Variable("u",false,range);
		Variable v=new Variable("v",false,range);
		
		//CsmDataReadStream  cdrs=new CsmDataReadStream(csd);
		//cdrs.readData(Type.CUBIC_P,u,v);
		//cdrs.closeFile();
		
		EliassenModelInCC tp=new EliassenModelInCC(csm);
		DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
		
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		Variable[] vel=ct.reprojectToCylindrical(u,v);
		
		Variable ut=vel[0];
		//Variable vr=vel[1];
		
		vel=ct.reprojectToCylindrical(u,v);
		Variable utL=vel[0]; utL.setName("utL");
		Variable vrL=vel[1]; vrL.setName("vrL");
		
		tp.cStormRelativeAziRadVelocity(utL,vrL);
		
		Variable aam1 =dm.cAbsoluteAngularMomentumByJohnson(ut);			aam1.setName("aam1");
		Variable aam2 =dm.cAbsoluteAngularMomentumByJohnsonWithSAA(ut);		aam2.setName("aam2");
		Variable aam3 =dm.cAbsoluteAngularMomentum(ut);						aam3.setName("aam3");
		Variable aam1L=dm.cAbsoluteAngularMomentumByJohnson(utL);			aam1L.setName("aam1L");
		Variable aam2L=dm.cAbsoluteAngularMomentumByJohnsonWithSAA(utL);	aam2L.setName("aam2L");
		Variable aam3L=dm.cAbsoluteAngularMomentum(utL);					aam3L.setName("aam3L");
		
		DataWrite dw=DataIOFactory.getDataWrite(csd.getCtlDescriptor(),"D:/Data/Validate/AAMDefinition/PAMcylind.dat");
		dw.writeData(csd.getCtlDescriptor(),aam1,aam2,aam3,aam1L,aam2L,aam3L);	dw.closeFile();
		
		Variable aam1m = aam1.anomalizeX(); aam1m.setName("aam1" );
		Variable aam2m = aam2.anomalizeX(); aam2m.setName("aam2" );
		Variable aam3m = aam3.anomalizeX(); aam3m.setName("aam3" );
		Variable aam1Lm=aam1L.anomalizeX();aam1Lm.setName("aam1L");
		Variable aam2Lm=aam2L.anomalizeX();aam2Lm.setName("aam2L");
		Variable aam3Lm=aam3L.anomalizeX();aam3Lm.setName("aam3L");
		
		Variable utm = ut.anomalizeX();
		Variable utLm=utL.anomalizeX();
		
		Variable aam4 =dm.cMeanAbsoluteAngularMomentumByJohnson(utm);	aam4.setName("aam4");
		Variable aam5 =dm.cMeanAbsoluteAngularMomentum(utm);			aam5.setName("aam5");
		Variable aam4L=dm.cMeanAbsoluteAngularMomentumByJohnson(utLm);	aam4L.setName("aam4L");
		Variable aam5L=dm.cMeanAbsoluteAngularMomentum(utLm);			aam5L.setName("aam5L");
		
		dw=DataIOFactory.getDataWrite(csd.getCtlDescriptor(),"D:/Data/Validate/AAMDefinition/PAM.dat");
		dw.writeData(csd.getCtlDescriptor(),utm,utLm,aam1m,aam2m,aam3m,aam4,aam5,aam1Lm,aam2Lm,aam3Lm,aam4L,aam5L);
		dw.closeFile();
	}
}
