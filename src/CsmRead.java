//
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class CsmRead{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/Validate/CsmRead/Haima.csm");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] vs=df.getVariables(new Range("t(3,3)",dd),"u","v");
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"D:/Data/Validate/CsmRead/cylindCP.dat");
		dw.writeData(dd,vs);
		dw.closeFile();
	}
}
