//
import miniufo.application.GeoFluidApplication.BoundaryCondition;
import miniufo.application.advanced.EllipticEqSORSolver;
import miniufo.application.advanced.EllipticEqSORSolver.DimCombination;
import miniufo.application.basic.DynamicMethodsInSC;
import miniufo.application.basic.SphericalHarmonicExpansion;
import miniufo.application.diagnosticModel.HGTInSC;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.DataWrite;


//
public class HGTModel{
	//
	private static final int   max=20000;
	private static final float tol=1e-4f;
	
	private static final String path="/lustre/home/qianyk/Data/";
	
	//
	public static void main(String[] args){
		ConcurrentUtil.initDefaultExecutor(9);
		
		//globalComputation(null,"global");
		regionalComputation("lon(90,240);lat(0,70)","WP");
		
		ConcurrentUtil.shutdown();
	}
	
	static void globalComputation(String hrange,String prefix){
		String range="";
		
		if(hrange!=null) range=";"+hrange; 
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"Haima.ctl");
		DataDescriptor ctl=df.getDataDescriptor();
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(ctl);
		HGTInSC gh=new HGTInSC(ssm);
		DynamicMethodsInSC gdm=new DynamicMethodsInSC(ssm);
		SphericalHarmonicExpansion sh=new SphericalHarmonicExpansion(ssm); 	sh.setM(180);
		EllipticEqSORSolver ees=new EllipticEqSORSolver(ssm);
		
		Variable[] vs=df.getVariables(new Range("z(35,37);t(1,3)"+range,ctl),"h","u","v","w");
		
		Variable h=vs[0].multiplyEq(9.8f);
		Variable u=vs[1];
		Variable v=vs[2];
		Variable w=vs[3];
		Variable div=gdm.c2DDivergence(u,v);
		Variable vor=gdm.c2DVorticity(u,v);
		
		float undef=h.getUndef();
		
		Variable f1 =gh.cTerm1(div);	f1.setUndef(undef);
		Variable f2 =gh.cTerm2(u);		f2.setUndef(undef);
		Variable f3 =gh.cTerm3(u,v);	f3.setUndef(undef);
		Variable f4 =gh.cTerm4(u,v);	f4.setUndef(undef);
		Variable f5 =gh.cTerm5(v);		f5.setUndef(undef);
		Variable f6 =gh.cTerm6(w,div);	f6.setUndef(undef);
		Variable f7 =gh.cTerm7(u,w);	f7.setUndef(undef);
		Variable f8 =gh.cTerm8(v,w);	f8.setUndef(undef);
		Variable f9 =gh.cTerm9(vor);	f9.setUndef(undef);
		Variable f10=gh.cTerm10(u);		f10.setUndef(undef);
		Variable f11=gh.cTerm11(u,v);	f11.setUndef(undef);
		Variable f12=gh.cTerm12(u);		f12.setUndef(undef);
		Variable ff=f1.plus(f2).plusEq(f3).plusEq(f4).plusEq(f5).
		plusEq(f6).plusEq(f7).plusEq(f8).plusEq(f9).plusEq(f10).plusEq(f11).plusEq(f12);
		ff.setName("ff");
		
		Variable h1 =(Variable)h.copy();	h1.setName("h1");	h1.setUndef(undef);		h1.setValue(0);
		Variable h2 =(Variable)h.copy();	h2.setName("h2");	h2.setUndef(undef);		h2.setValue(0);
		Variable h3 =(Variable)h.copy();	h3.setName("h3");	h3.setUndef(undef);		h3.setValue(0);
		Variable h4 =(Variable)h.copy();	h4.setName("h4");	h4.setUndef(undef);		h4.setValue(0);
		Variable h5 =(Variable)h.copy();	h5.setName("h5");	h5.setUndef(undef);		h5.setValue(0);
		Variable h6 =(Variable)h.copy();	h6.setName("h6");	h6.setUndef(undef);		h6.setValue(0);
		Variable h7 =(Variable)h.copy();	h7.setName("h7");	h7.setUndef(undef);		h7.setValue(0);
		Variable h8 =(Variable)h.copy();	h8.setName("h8");	h8.setUndef(undef);		h8.setValue(0);
		Variable h9 =(Variable)h.copy();	h9.setName("h9");	h9.setUndef(undef);		h9.setValue(0);
		Variable h10=(Variable)h.copy();	h10.setName("h10");	h10.setUndef(undef);	h10.setValue(0);
		Variable h11=(Variable)h.copy();	h11.setName("h11");	h11.setUndef(undef);	h11.setValue(0);
		Variable h12=(Variable)h.copy();	h12.setName("h12");	h12.setUndef(undef);	h12.setValue(0);
		Variable hb =(Variable)h.copy();	hb.setName("hb");	hb.setUndef(undef);		hb.setInner(0);	
		Variable hh =(Variable)h.copy();	hh.setName("hh");	hh.setUndef(undef);		hh.setInner(0);	
		
		ees.setDimCombination(DimCombination.XY);
		ees.setBCofX(BoundaryCondition.Periodic);
		ees.setBCofY(BoundaryCondition.Fixed);
		ees.setTolerance(tol);
		ees.setMaxLoopCount(max);
		
		ees.solve(h1 ,f1 );	ees.solve(h8 ,f8 );
		ees.solve(h2 ,f2 );	ees.solve(h9 ,f9 );
		ees.solve(h3 ,f3 );	ees.solve(h10,f10);
		ees.solve(h4 ,f4 );	ees.solve(h11,f11);
		ees.solve(h5 ,f5 );	ees.solve(h12,f12);
		ees.solve(h6 ,f6 );	ees.solve(hh ,ff);
		ees.solve(h7 ,f7 );	ees.solve(hb ,null);
		
		Variable h1sh =sh.solvePoissonEquation(f1 );	h1sh.setName("h1sh");	h1sh.setUndef(undef);
		Variable h2sh =sh.solvePoissonEquation(f2 );	h2sh.setName("h2sh");	h2sh.setUndef(undef);
		Variable h3sh =sh.solvePoissonEquation(f3 );	h3sh.setName("h3sh");	h3sh.setUndef(undef);
		Variable h4sh =sh.solvePoissonEquation(f4 );	h4sh.setName("h4sh");	h4sh.setUndef(undef);
		Variable h5sh =sh.solvePoissonEquation(f5 );	h5sh.setName("h5sh");	h5sh.setUndef(undef);
		Variable h6sh =sh.solvePoissonEquation(f6 );	h6sh.setName("h6sh");	h6sh.setUndef(undef);
		Variable h7sh =sh.solvePoissonEquation(f7 );	h7sh.setName("h7sh");	h7sh.setUndef(undef);
		Variable h8sh =sh.solvePoissonEquation(f8 );	h8sh.setName("h8sh");	h8sh.setUndef(undef);
		Variable h9sh =sh.solvePoissonEquation(f9 );	h9sh.setName("h9sh");	h9sh.setUndef(undef);
		Variable h10sh=sh.solvePoissonEquation(f10);	h10sh.setName("h10sh");	h10sh.setUndef(undef);
		Variable h11sh=sh.solvePoissonEquation(f11);	h11sh.setName("h11sh");	h11sh.setUndef(undef);
		Variable h12sh=sh.solvePoissonEquation(f12);	h12sh.setName("h12sh");	h12sh.setUndef(undef);
		Variable hhsh =sh.solvePoissonEquation(ff );	hhsh.setName("hhsh");	hhsh.setUndef(undef);
		
		DataWrite dw=new miniufo.io.CtlDataWriteStream(path+prefix+"HGT.dat");
		dw.writeData(ctl,h,u,v,w,div,vor,hh,hb,ff,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,
			h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,
			h1sh,h2sh,h3sh,h4sh,h5sh,h6sh,h7sh,h8sh,h9sh,h10sh,h11sh,h12sh,hhsh
		);
		dw.closeFile();
	}
	
	static void regionalComputation(String hrange,String prefix){
		String range="";
		
		if(hrange!=null) range=";"+hrange; 
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"Haima.ctl");
		DataDescriptor ctl=df.getDataDescriptor();
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(ctl);
		HGTInSC gh=new HGTInSC(ssm);
		DynamicMethodsInSC gdm=new DynamicMethodsInSC(ssm);
		EllipticEqSORSolver ees=new EllipticEqSORSolver(ssm);
		
		Variable[] vs=df.getVariables(new Range("z(35,37);t(1,3)"+range,ctl),"h","u","v","w");
		
		Variable h=vs[0].multiplyEq(9.8f);
		Variable u=vs[1];
		Variable v=vs[2];
		Variable w=vs[3];
		Variable div=gdm.c2DDivergence(u,v);
		Variable vor=gdm.c2DVorticity(u,v);
		
		float undef=h.getUndef();
		
		Variable f1 =gh.cTerm1(div);	f1.setUndef(undef);
		Variable f2 =gh.cTerm2(u);		f2.setUndef(undef);
		Variable f3 =gh.cTerm3(u,v);	f3.setUndef(undef);
		Variable f4 =gh.cTerm4(u,v);	f4.setUndef(undef);
		Variable f5 =gh.cTerm5(v);		f5.setUndef(undef);
		Variable f6 =gh.cTerm6(w,div);	f6.setUndef(undef);
		Variable f7 =gh.cTerm7(u,w);	f7.setUndef(undef);
		Variable f8 =gh.cTerm8(v,w);	f8.setUndef(undef);
		Variable f9 =gh.cTerm9(vor);	f9.setUndef(undef);
		Variable f10=gh.cTerm10(u);		f10.setUndef(undef);
		Variable f11=gh.cTerm11(u,v);	f11.setUndef(undef);
		Variable f12=gh.cTerm12(u);		f12.setUndef(undef);
		Variable ff=f1.plus(f2).plusEq(f3).plusEq(f4).plusEq(f5).
		plusEq(f6).plusEq(f7).plusEq(f8).plusEq(f9).plusEq(f10).plusEq(f11).plusEq(f12);
		ff.setName("ff");
		
		Variable h1 =(Variable)h.copy();	h1.setName("h1");	h1.setUndef(undef);		h1.setValue(0);
		Variable h2 =(Variable)h.copy();	h2.setName("h2");	h2.setUndef(undef);		h2.setValue(0);
		Variable h3 =(Variable)h.copy();	h3.setName("h3");	h3.setUndef(undef);		h3.setValue(0);
		Variable h4 =(Variable)h.copy();	h4.setName("h4");	h4.setUndef(undef);		h4.setValue(0);
		Variable h5 =(Variable)h.copy();	h5.setName("h5");	h5.setUndef(undef);		h5.setValue(0);
		Variable h6 =(Variable)h.copy();	h6.setName("h6");	h6.setUndef(undef);		h6.setValue(0);
		Variable h7 =(Variable)h.copy();	h7.setName("h7");	h7.setUndef(undef);		h7.setValue(0);
		Variable h8 =(Variable)h.copy();	h8.setName("h8");	h8.setUndef(undef);		h8.setValue(0);
		Variable h9 =(Variable)h.copy();	h9.setName("h9");	h9.setUndef(undef);		h9.setValue(0);
		Variable h10=(Variable)h.copy();	h10.setName("h10");	h10.setUndef(undef);	h10.setValue(0);
		Variable h11=(Variable)h.copy();	h11.setName("h11");	h11.setUndef(undef);	h11.setValue(0);
		Variable h12=(Variable)h.copy();	h12.setName("h12");	h12.setUndef(undef);	h12.setValue(0);
		Variable hb =(Variable)h.copy();	hb.setName("hb");	hb.setUndef(undef);		hb.setInner(0);	
		Variable hh =(Variable)h.copy();	hh.setName("hh");	hh.setUndef(undef);		hh.setInner(0);	
		
		ees.setDimCombination(DimCombination.XY);
		ees.setBCofX(BoundaryCondition.Periodic);
		ees.setBCofY(BoundaryCondition.Fixed);
		ees.setTolerance(tol);
		ees.setMaxLoopCount(max);
		
		ees.solve(h1 ,f1 );	ees.solve(h8 ,f8 );
		ees.solve(h2 ,f2 );	ees.solve(h9 ,f9 );
		ees.solve(h3 ,f3 );	ees.solve(h10,f10);
		ees.solve(h4 ,f4 );	ees.solve(h11,f11);
		ees.solve(h5 ,f5 );	ees.solve(h12,f12);
		ees.solve(h6 ,f6 );	ees.solve(hh ,ff);
		ees.solve(h7 ,f7 );	ees.solve(hb ,null);
		
		DataWrite dw=new miniufo.io.CtlDataWriteStream(path+prefix+"HGTEB.dat");
		dw.writeData(ctl,h,u,v,w,div,vor,hh,hb,ff,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,
			h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12
		);
		dw.closeFile();
	}
}
