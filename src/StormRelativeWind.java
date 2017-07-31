//
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.CsmDataReadStream;
import miniufo.io.CsmDataWriteStream;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class StormRelativeWind{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.csm");
		CsmDescriptor csd=(CsmDescriptor)df.getDataDescriptor();
		CylindricalSpatialModel csm=(CylindricalSpatialModel)df.getSpatialModel();
		SphericalSpatialModel ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		
		Range range=new Range("",csd);
		
		Variable uo=new Variable("u",false,range);
		Variable vo=new Variable("v",false,range);
		
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd);
		cdrs.readData(Type.CUBIC_P,uo,vo); cdrs.closeFile();
		
		EliassenModelInCC tp=new EliassenModelInCC(csm);
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		Variable uL=uo.copy();
		Variable vL=vo.copy();
		
		// method 1: storm-relative and then reproject
		tp.cStormRelativeLatLonVelocity(uL,vL);
		Variable[] vel=ct.reprojectToCylindrical(uL,vL);
		Variable utL1=vel[0];
		Variable vrL1=vel[1];
		
		// method 2: reproject and then storm-relative
		vel=ct.reprojectToCylindrical(uo,vo);
		Variable ut=vel[0].copy();
		Variable vr=vel[1].copy();
		Variable utL2=vel[0]; Variable utm=vel[0].copy().anomalizeX();
		Variable vrL2=vel[1]; Variable vrm=vel[1].copy().anomalizeX();
		tp.cStormRelativeAziRadVelocity(utL2,vrL2);
		
		ut.setName("ut");
		vr.setName("vr");
		utL1.setName("utL1");
		vrL1.setName("vrL1");
		utL2.setName("utL2");
		vrL2.setName("vrL2");
		
		vel=ct.reprojectToLatLon(utL1,vrL1);
		Variable uL1=vel[0]; uL1.setName("uL1");
		Variable vL1=vel[1]; vL1.setName("vL1");
		
		vel=ct.reprojectToLatLon(utL2,vrL2);
		Variable uL2=vel[0]; uL2.setName("uL2");
		Variable vL2=vel[1]; vL2.setName("vL2");
		
		CsmDataWriteStream cdws=new CsmDataWriteStream("d:/SRCylin.dat");
		cdws.writeData(csd,uo,vo,ut,vr,utL1,vrL1,utL2,vrL2,uL1,vL1,uL2,vL2);
		
		Variable utLm1=utL1.anomalizeX();
		Variable utLm2=utL2.anomalizeX();
		Variable vrLm1=vrL1.anomalizeX();
		Variable vrLm2=vrL2.anomalizeX();
		
		DataWrite dw=DataIOFactory.getDataWrite(csd.getCtlDescriptor(),"D:/SR.dat");
		dw.writeData(csd.getCtlDescriptor(),utm,vrm,utLm1,vrLm1,utLm2,vrLm2);	dw.closeFile();
	}
}
