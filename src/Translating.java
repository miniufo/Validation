//
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.CsmDataWriteStream;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class Translating{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.csm");
		CsmDescriptor csd=(CsmDescriptor)df.getDataDescriptor();
		CylindricalSpatialModel csm=new CylindricalSpatialModel(csd);
		SphericalSpatialModel ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		
		EliassenModelInCC tp=new EliassenModelInCC(csm);
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		Variable[] whole1=null;//tp.whole1;
		Variable[] whole2=null;//tp.whole2;
		Variable[] whole3=null;//tp.whole3;
		
		Variable[] vel=null;
		
		vel=ct.reprojectToCylindrical(whole1[0],whole1[1]);
		Variable ut1=vel[0],vr1=vel[1];
		
		vel=ct.reprojectToCylindrical(whole2[0],whole2[1]);
		Variable ut2=vel[0],vr2=vel[1];
		
		vel=ct.reprojectToCylindrical(whole3[0],whole3[1]);
		Variable ut3=vel[0],vr3=vel[1];
		
		whole1[0].setName("u1"); whole1[1].setName("v1");
		whole2[0].setName("u2"); whole2[1].setName("v2");
		whole3[0].setName("u3"); whole3[1].setName("v3");
		
		ut1.setName("ut1"); vr1.setName("vr1");
		ut2.setName("ut2"); vr2.setName("vr2");
		ut3.setName("ut3"); vr3.setName("vr3");
		
		CsmDataWriteStream cdws=new CsmDataWriteStream("D:/Data/Validate/StormRelative/TransCylin.dat");
		cdws.writeData(csd,whole1[0],whole1[1],whole2[0],whole2[1],whole3[0],whole3[1],ut1,ut2,ut3,vr1,vr2,vr3);
		
		Variable utm1=ut1.anomalizeX(); Variable um1=whole1[0].anomalizeX();
		Variable vrm1=vr1.anomalizeX(); Variable vm1=whole1[1].anomalizeX();
		Variable utm2=ut2.anomalizeX(); Variable um2=whole2[0].anomalizeX();
		Variable vrm2=vr2.anomalizeX(); Variable vm2=whole2[1].anomalizeX();
		Variable utm3=ut3.anomalizeX(); Variable um3=whole3[0].anomalizeX();
		Variable vrm3=vr3.anomalizeX(); Variable vm3=whole3[1].anomalizeX();
		
		DataWrite dw=DataIOFactory.getDataWrite(csd.getCtlDescriptor(),"D:/Data/Validate/StormRelative/Trans.dat");
		dw.writeData(csd.getCtlDescriptor(),utm1,vrm1,utm2,vrm2,utm3,vrm3,um1,um2,um3,vm1,vm2,vm3);	dw.closeFile();
	}
}
