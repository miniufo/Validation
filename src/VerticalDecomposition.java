//
import miniufo.application.basic.ThermoDynamicMethodsInSC;
import miniufo.application.statisticsModel.EOFApplication;
import miniufo.basic.ArrayUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class VerticalDecomposition{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("/lustre/home/qianyk/Data/Haima.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		ThermoDynamicMethodsInSC tm=new ThermoDynamicMethodsInSC(ssm);
		
		Variable[] vs=df.getVariables(new Range("lon(122,122);lat(10,40);t(1,3)",dd),"t","u");
		
		Variable T=vs[0];
		Variable u=vs[1];
		Variable S=tm.cStaticStabilityArgByT(T);
		
		Variable[] modes=EOFApplication.verticalDecompose(S,u,37,dd.getDZDef()[0]);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"/lustre/home/qianyk/Data/VMD.dat");
		dw.writeData(dd,ArrayUtil.concatAll(Variable.class,new Variable[]{T,u,S},modes));
		dw.closeFile();
	}
}
