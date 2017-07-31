//
import miniufo.application.statisticsModel.FilterMethods;
import miniufo.application.statisticsModel.WaveNumberFrequencyAnalysis;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.mathsphysics.WindowFunction;


//
public class WaveNoFreq{
	
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/Validate/WaveNoFreq/hgt.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable hgt=df.getVariables(new Range("t(1,14600)",dd),"hgt")[0];
		
		float mean=hgt.averageAlong(Dimension.T).averageAlong(Dimension.X).getData()[0][0][0][0];
		
		System.out.println(mean);
		
		hgt.minusEq(mean);
		
		FilterMethods.FourierFilter(hgt,Dimension.T,1460);
		FilterMethods.FourierFilter(hgt,Dimension.T,730);
		FilterMethods.FourierFilter(hgt,Dimension.T,365);
		FilterMethods.removeLinearTrend(hgt);
		
		Variable h=new Variable("hgt",false,new Range(14000,1,1,hgt.getXCount()));
		
		for(int i=0;i<hgt.getXCount();i++)
		System.arraycopy(hgt.getData()[0][0][i],0,h.getData()[0][0][i],0,h.getTCount());
		
		WaveNumberFrequencyAnalysis wnfa=new WaveNumberFrequencyAnalysis(h);
		wnfa.windowingT(WindowFunction.tukey(h.getTCount(),0.05f));
		wnfa.transform();
		
		Variable mod =wnfa.getMode();
		Variable real=wnfa.getRealPart();
		Variable imag=wnfa.getImagPart();
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/Data/Validate/WaveNoFreq/wave.dat");
		dw.writeData(dd,h,mod,real,imag);	dw.closeFile();
	}
}
