//
import miniufo.application.basic.VelocityFieldInSC;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class VeloctiyField{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/Validate/VectorField/uv.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		VelocityFieldInSC vf=new VelocityFieldInSC(ssm);
		
		Variable[] vs=df.getVariables(new Range("",dd),"u","v");
		
		Variable[] uv=vf.cRotationalVelocityByEndlich(vs[0],vs[1]);
		
		Variable sf=vf.cStreamFunctionByEndlich(uv[0],uv[1]);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/Data/Validate/VectorField/sf.dat");
		dw.writeData(dd,uv[0],uv[1],vs[0],vs[1],sf);
		dw.closeFile();
	}
}
