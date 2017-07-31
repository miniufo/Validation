//
import miniufo.application.basic.DynamicMethodsInSC;
import miniufo.test.application.basic.GlobalLaplaceEquationInSC;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class StrainRate{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		DynamicMethodsInSC dm=new DynamicMethodsInSC(ssm);
		GlobalLaplaceEquationInSC le=new GlobalLaplaceEquationInSC(ssm);
		
		Variable[] vs=df.getVariables(new Range("t(1,1);lev(500,500)",dd),"u","v");
		
		Variable tstr=dm.c2DTensionStrain(vs[0],vs[1]);
		Variable sstr=dm.c2DShearingStrain(vs[0],vs[1]);
		Variable vor =dm.c2DVorticity(vs[0],vs[1]);
		Variable div =dm.c2DDivergence(vs[0],vs[1]);
		
		Variable sf=new Variable("sf",vor);	le.solve(5000,sf,vor);
		Variable pf=new Variable("pf",vor);	le.solve(5000,pf,div);
		Variable sf1=new Variable("sf1",vor);	le.solve(5000,sf1,tstr);
		Variable sf2=new Variable("sf2",vor);	le.solve(5000,sf2,sstr);
		
		Variable[] grdsf=dm.c2DGradient(sf);	grdsf[0].setName("sfx");	grdsf[1].setName("sfy");
		Variable[] grdpf=dm.c2DGradient(pf);	grdpf[0].setName("pfx");	grdpf[1].setName("pfy");
		Variable[] grdsf1=dm.c2DGradient(sf1);	grdsf1[0].setName("sf1x");	grdsf1[1].setName("sf1y");
		Variable[] grdsf2=dm.c2DGradient(sf2);	grdsf2[0].setName("sf2x");	grdsf2[1].setName("sf2y");
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/Data/Validate/StrainRate/strain.dat");
		dw.writeData(dd,vs[0],vs[1],tstr,sstr,vor,div,sf,pf,sf1,sf2,grdsf[0],grdsf[1],grdpf[0],grdpf[1],grdsf1[0],grdsf1[1],grdsf2[0],grdsf2[1]);
		dw.closeFile();
	}
}
