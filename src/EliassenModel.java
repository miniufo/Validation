//
import miniufo.application.basic.DynamicMethodsInSC;
import miniufo.application.basic.ThermoDynamicMethodsInSC;
import miniufo.application.diagnosticModel.EliassenModelInSC;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.util.DataInterpolation;


public final class EliassenModel{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/ERAInterim/Eliassen/Data.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		
		EliassenModelInSC emd=new EliassenModelInSC(ssm);
		ThermoDynamicMethodsInSC tmd=new ThermoDynamicMethodsInSC(ssm);
		ThermoDynamicMethodsInSC tdmd=new ThermoDynamicMethodsInSC(ssm);
		DynamicMethodsInSC dm=new DynamicMethodsInSC(ssm);
		
		Variable[] vs=df.getVariables(new Range("t(1,1)",dd),"z","t","q","w","u","v");
		
		Variable z=vs[0];
		Variable t=vs[1];
		Variable q=vs[2];
		Variable w=vs[3];
		Variable u=vs[4];
		Variable v=vs[5];
		
		Variable s =tmd.cStaticStabilityArgByPT(tmd.cPotentialTemperature(t));
		Variable g =emd.cAbsoluteAngularMomentum(u);
		Variable Q =tdmd.cLargeScaleLatentHeating(q,u,v,w,t);
		Variable th=tmd.cPotentialTemperature(t);
		
		Variable zm = z.anomalizeX();
		Variable tm = t.anomalizeX();
		Variable qm = q.anomalizeX();
		Variable wm = w.anomalizeX();
		Variable um = u.anomalizeX();
		Variable vm = v.anomalizeX();
		Variable sm = s.anomalizeX();
		Variable gm = g.anomalizeX();
		Variable Qm = Q.anomalizeX();
		Variable thm=th.anomalizeX();
		
		Variable gava= g.multiply(v).anomalizeX(); gava.setName("gava");
		Variable gawa= g.multiply(w).anomalizeX(); gawa.setName("gawa");
		Variable tava=th.multiply(v).anomalizeX(); tava.setName("tava");
		Variable tawa=th.multiply(w).anomalizeX(); tawa.setName("tawa");
		
		Variable sfe=emd.cEddyInducedStreamfunction(tava,thm);
		Variable[] vedy=emd.cEddyInducedVelocity(sfe);
		Variable[] EP=emd.cEPVector(tava,gava,gawa,gm,thm,0.8f);
		Variable div=dm.cYZDivergence(EP[0],EP[1]);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"D:/Data/ERAInterim/Eliassen/Hadley.dat");
		dw.writeData(dd,zm,tm,qm,wm,um,vm,sm,gm,Qm,thm,gava,gawa,tava,tawa,sfe,vedy[0],vedy[1],EP[0],EP[1],EP[2],EP[3],div);
		dw.closeFile();
	}
	
	static void vInterp(){
		DiagnosisFactory df=DiagnosisFactory.parseFile("D:/Data/ERAInterim/Eliassen/DataOrig.ctl");
		DataInterpolation di=new DataInterpolation(df.getDataDescriptor());
		di.verticalInterp("d:/Data/ERAInterim/Eliassen/Data.dat",39);
	}
}
