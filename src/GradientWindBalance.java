//
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.CsmDataReadStream;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class GradientWindBalance{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.csm");
		CsmDescriptor csd=(CsmDescriptor)df.getDataDescriptor();
		CylindricalSpatialModel csm=(CylindricalSpatialModel)df.getSpatialModel();
		SphericalSpatialModel ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		
		Range range=new Range("",csd);
		
		Variable u=new Variable("u",false,range);
		Variable v=new Variable("v",false,range);
		Variable h=new Variable("h",false,range);
		
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd);
		cdrs.readData(Type.CUBIC_P,u,v,h); h.multiplyEq(9.8f);
		cdrs.closeFile();
		
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
		
		//tp.cStormRelativeCylindricalVelocity(utL,vrL);
		
		Variable utm = ut.anomalizeX();
		Variable utLm=utL.anomalizeX();
		Variable hm  =  h.anomalizeX();
		
		Variable gw1=dm.cMeanGradientWindByJohnson(hm);	gw1.setName("gw1");
		Variable gw2=dm.cMeanGradientWindByCurvedGWB(hm);	gw2.setName("gw2");
		
		DataWrite dw=DataIOFactory.getDataWrite(csd.getCtlDescriptor(),"D:/Data/Validate/GradientWindBalance/GWB.dat");
		dw.writeData(csd.getCtlDescriptor(),utm,utLm,gw1,gw2);	dw.closeFile();
	}
	
	static void cAAM(){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.csm");
		CsmDescriptor csd=(CsmDescriptor)df.getDataDescriptor();
		CylindricalSpatialModel csm=(CylindricalSpatialModel)df.getSpatialModel();
		SphericalSpatialModel ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		
		Range range=new Range("",csd);
		
		Variable u=new Variable("u",false,range);
		Variable v=new Variable("v",false,range);
		Variable h=new Variable("h",false,range);
		
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd);
		cdrs.readData(Type.CUBIC_P,u,v,h); h.multiplyEq(9.8f);
		cdrs.closeFile();
		
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
		
		Variable aam1 =dm.cAbsoluteAngularMomentumByJohnson(ut);
		Variable aam2 =dm.cAbsoluteAngularMomentumByJohnsonWithSAA(ut);
		Variable aam3 =dm.cAbsoluteAngularMomentum(ut);
		Variable aam1L=dm.cAbsoluteAngularMomentumByJohnson(utL);
		Variable aam2L=dm.cAbsoluteAngularMomentumByJohnsonWithSAA(utL);
		Variable aam3L=dm.cAbsoluteAngularMomentum(utL);
		
		Variable aam1m = aam1.anomalizeX(); aam1m.setName("aam1" );
		Variable aam2m = aam2.anomalizeX(); aam2m.setName("aam2" );
		Variable aam3m = aam3.anomalizeX(); aam3m.setName("aam3" );
		Variable aam1Lm=aam1L.anomalizeX();aam1Lm.setName("aam1L");
		Variable aam2Lm=aam2L.anomalizeX();aam2Lm.setName("aam2L");
		Variable aam3Lm=aam3L.anomalizeX();aam3Lm.setName("aam3L");
		
		Variable utm = ut.anomalizeX();
		Variable utLm=utL.anomalizeX();
		
		DataWrite dw=DataIOFactory.getDataWrite(csd.getCtlDescriptor(),"D:/Data/Validate/GradientWindBalance/AAM.dat");
		dw.writeData(csd.getCtlDescriptor(),utm,utLm,aam1m,aam2m,aam3m,aam1Lm,aam2Lm,aam3Lm);	dw.closeFile();
	}
}
