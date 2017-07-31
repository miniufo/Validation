//
import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.contour.ContourSphericalSpatialModel;
import miniufo.application.contour.KeffInCTS;
import miniufo.application.contour.KeffInSC;
import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.CtsDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.util.DataInterpolation;

//
public final class KeffTest{
	//
	private static final String path="D:/Data/Validate/Keff/";
	
	//
	public static void main(String[] args){
		//generateTracers(); System.exit(0);
		
		//cCartesianKeffs("Linear",200);
		//cCartesianKeffs("Linear2",400);
		//cCartesianKeffs("Linear4",800);
		//cCartesianKeffs("Sin",200);
		//cCartesianKeffs("Sin2",400);
		//cCartesianKeffs("Sin4",800);
		
		//cSphericalKeffs("Linear",181);
		//cSphericalKeffs("Linear2",361);
		//cSphericalKeffs("Linear4",721);
		//cSphericalKeffs("Sin",181);
		//cSphericalKeffs("Sin2",361);
		//cSphericalKeffs("Sin4",721);
		
		//interpPVData(2);
		cPVKeffs("PV" ,241);
		//cPVKeffs("PV2",481);
		//cPVKeffs("PV4",961);
		//cPVKeffs("PV8",1921);
	}
	
	
	static void interpPVData(int res){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"PVTrue/PV.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		DataInterpolation di=new DataInterpolation(dd);
		di.horizontalInterp(path+"PVTrue/PV"+res+".dat",Type.PERIODIC_LINEAR,Type.LINEAR,(dd.getYCount()-1)*res+1,dd.getXCount()*res);
	}
	
	static void cPVKeffs(String fname,int interpY){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"PVTrue/"+fname+".ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] v0=cPVKeff(df, 41,interpY);
		Variable[] v1=cPVKeff(df, 81,interpY);
		Variable[] v2=cPVKeff(df,161,interpY);
		Variable[] v3=cPVKeff(df,321,interpY);
		
		for(Variable v:v0) v.setName(v.getName()+"c41" );
		for(Variable v:v1) v.setName(v.getName()+"c81" );
		for(Variable v:v2) v.setName(v.getName()+"c161");
		for(Variable v:v3) v.setName(v.getName()+"c321");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"PVTrue/"+fname+"Keffs.dat");
		cdws.writeData(dd,ArrayUtil.concatAll(Variable.class,v0,v1,v2,v3));
		cdws.closeFile();
	}
	
	static Variable[] cPVKeff(DiagnosisFactory df,int qnum,int interpY){
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable pv=df.getVariables(new Range("t(1,1)",dd),"pv")[0].multiplyEq(1e6f);
		
		ContourSphericalSpatialModel cc=new ContourSphericalSpatialModel(dd);
		cc.initContourByTracer(pv,qnum,1,true);System.out.println(cc.getContours()[0][0]);
		
		KeffInSC keffSC=new KeffInSC(cc);
		
		Variable area   =cc.getAreasBoundedByContour();
		Variable intGrd2=cc.integrateWithinContour(cc.getSquaredTracerGradient());
		
		Variable aveGrd2AlgC=keffSC.cGradientWRTArea(intGrd2);
		Variable qGrdA  =keffSC.cGradientWRTArea();
		Variable dqdye  =keffSC.cDqDye();
		Variable Le2    =keffSC.cEquivalentLengthSquare(aveGrd2AlgC,qGrdA);
		Variable Lmin2  =keffSC.cMinimumLengthSquare();
		Variable nkeff  =keffSC.cNormalizedKeff(aveGrd2AlgC,dqdye);
		
		Type t=Type.LINEAR;
		
		area       =cc.interpolatedToYs(area       ,interpY,t);
		qGrdA      =cc.interpolatedToYs(qGrdA      ,interpY,t);
		intGrd2    =cc.interpolatedToYs(intGrd2    ,interpY,t);
		aveGrd2AlgC=cc.interpolatedToYs(aveGrd2AlgC,interpY,t);
		dqdye      =cc.interpolatedToYs(dqdye      ,interpY,t);
		Le2        =cc.interpolatedToYs(Le2        ,interpY,t);
		Lmin2      =cc.interpolatedToYs(Lmin2      ,interpY,t);
		nkeff      =cc.interpolatedToYs(nkeff      ,interpY,t);
		
		return new Variable[]{area,qGrdA,intGrd2,aveGrd2AlgC,dqdye,nkeff,Le2,Lmin2};
	}
	
	
	static void cCartesianKeffs(String fname,int interpY){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"Cartesian/"+fname+".cts");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] v0=cCartesianKeff(df, 11,interpY);
		Variable[] v1=cCartesianKeff(df, 31,interpY);
		Variable[] v2=cCartesianKeff(df, 61,interpY);
		Variable[] v3=cCartesianKeff(df,101,interpY);
		
		for(Variable v:v0) v.setName(v.getName()+"c11" );
		for(Variable v:v1) v.setName(v.getName()+"c31" );
		for(Variable v:v2) v.setName(v.getName()+"c61" );
		for(Variable v:v3) v.setName(v.getName()+"c101");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Cartesian/"+fname+"Keffs.dat");
		cdws.writeData(dd,ArrayUtil.concatAll(Variable.class,v0,v1,v2,v3));
		cdws.closeFile();
	}
	
	static Variable[] cCartesianKeff(DiagnosisFactory df,int qnum,int interpY){
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable pv=df.getVariables(new Range("t(1,1)",dd),"q")[0];
		
		ContourCartesianSpatialModel cc=new ContourCartesianSpatialModel((CtsDescriptor)dd);
		cc.initContourByTracer(pv,qnum,1,true);System.out.println(cc.getContours()[0][0]);
		
		KeffInCTS keffSC=new KeffInCTS(cc);
		
		Variable area   =cc.getAreasBoundedByContour();
		Variable intGrd2=cc.integrateWithinContour(cc.getSquaredTracerGradient());
		
		Variable aveGrd2AlgC=keffSC.cGradientWRTArea(intGrd2);
		Variable qGrdA  =keffSC.cGradientWRTArea();
		Variable dqdye  =keffSC.cDqDye();
		Variable Le2    =keffSC.cEquivalentLengthSquare(aveGrd2AlgC,qGrdA);
		Variable Lmin2  =keffSC.cMinimumLengthSquare();
		Variable nkeff  =keffSC.cNormalizedKeff(aveGrd2AlgC,dqdye);
		
		Type t=Type.LINEAR;
		
		area       =cc.interpolatedToYs(area       ,interpY,t);
		qGrdA      =cc.interpolatedToYs(qGrdA      ,interpY,t);
		intGrd2    =cc.interpolatedToYs(intGrd2    ,interpY,t);
		aveGrd2AlgC=cc.interpolatedToYs(aveGrd2AlgC,interpY,t);
		dqdye      =cc.interpolatedToYs(dqdye      ,interpY,t);
		Le2        =cc.interpolatedToYs(Le2        ,interpY,t);
		Lmin2      =cc.interpolatedToYs(Lmin2      ,interpY,t);
		nkeff      =cc.interpolatedToYs(nkeff      ,interpY,t);
		
		return new Variable[]{area,qGrdA,intGrd2,aveGrd2AlgC,dqdye,nkeff,Le2,Lmin2};
	}
	
	
	static void cSphericalKeffs(String fname,int interpY){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"Spherical/"+fname+".ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] v0=cSphericalKeff(df, 11,interpY);
		Variable[] v1=cSphericalKeff(df, 31,interpY);
		Variable[] v2=cSphericalKeff(df, 61,interpY);
		Variable[] v3=cSphericalKeff(df,101,interpY);
		
		for(Variable v:v0) v.setName(v.getName()+"c11" );
		for(Variable v:v1) v.setName(v.getName()+"c31" );
		for(Variable v:v2) v.setName(v.getName()+"c61" );
		for(Variable v:v3) v.setName(v.getName()+"c101");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Spherical/"+fname+"Keffs.dat");
		cdws.writeData(dd,ArrayUtil.concatAll(Variable.class,v0,v1,v2,v3));
		cdws.closeFile();
	}
	
	static Variable[] cSphericalKeff(DiagnosisFactory df,int qnum,int interpY){
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable pv=df.getVariables(new Range("t(1,1)",dd),"q")[0];
		
		ContourSphericalSpatialModel cc=new ContourSphericalSpatialModel(dd);
		cc.initContourByTracer(pv,qnum,1,true);System.out.println(cc.getContours()[0][0]);
		
		KeffInSC keffSC=new KeffInSC(cc);
		
		Variable area   =cc.getAreasBoundedByContour();
		Variable intGrd2=cc.integrateWithinContour(cc.getSquaredTracerGradient());
		
		Variable aveGrd2AlgC=keffSC.cGradientWRTArea(intGrd2);
		Variable qGrdA  =keffSC.cGradientWRTArea();
		Variable dqdye  =keffSC.cDqDye();
		Variable Le2    =keffSC.cEquivalentLengthSquare(aveGrd2AlgC,qGrdA);
		Variable Lmin2  =keffSC.cMinimumLengthSquare();
		Variable nkeff  =keffSC.cNormalizedKeff(aveGrd2AlgC,dqdye);
		
		Type t=Type.LINEAR;
		
		area       =cc.interpolatedToYs(area       ,interpY,t);
		qGrdA      =cc.interpolatedToYs(qGrdA      ,interpY,t);
		intGrd2    =cc.interpolatedToYs(intGrd2    ,interpY,t);
		aveGrd2AlgC=cc.interpolatedToYs(aveGrd2AlgC,interpY,t);
		dqdye      =cc.interpolatedToYs(dqdye      ,interpY,t);
		Le2        =cc.interpolatedToYs(Le2        ,interpY,t);
		Lmin2      =cc.interpolatedToYs(Lmin2      ,interpY,t);
		nkeff      =cc.interpolatedToYs(nkeff      ,interpY,t);
		
		return new Variable[]{area,qGrdA,intGrd2,aveGrd2AlgC,dqdye,nkeff,Le2,Lmin2};
	}
	
	
	static void generateTracers(){
		/**/
		genCartesianTracer(true ,1);
		genCartesianTracer(true ,2);
		genCartesianTracer(true ,4);
		genCartesianTracer(false,1);
		genCartesianTracer(false,2);
		genCartesianTracer(false,4);
		
		genSphericalTracer(true ,1);
		genSphericalTracer(true ,2);
		genSphericalTracer(true ,4);
		genSphericalTracer(false,1);
		genSphericalTracer(false,2);
		genSphericalTracer(false,4);
	}
	
	static void genSphericalTracer(boolean linear,int res){
		if(res<1) throw new IllegalArgumentException("res should be positive");
		
		int oy=181,ox=360;
		int y=(oy-1)*res+1,x=ox*res;
		
		Variable tr=new Variable("tr",new Range(1,1,y,x));
		
		float[][] tdata=tr.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			if(linear)
				tdata[j][i]=(float)(-1.0+2.0*j/(y-1.0)); // Linear
			else
				tdata[j][i]=(float)Math.sin(Math.PI/2.0*(-1.0+2.0*j/(y-1.0))); // Sin
		}
		
		String fname=linear?"Linear":"Sin";
		
		if(res!=1) fname+=res;
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Spherical/"+fname+".dat");
		cdws.writeData(tr); cdws.closeFile();
	}
	
	static void genCartesianTracer(boolean linear,int res){
		if(res<1) throw new IllegalArgumentException("res should be positive");
		
		int oy=200,ox=280;
		int y=(oy-1)*res+1,x=(ox-1)*res+1;
		
		Variable tr=new Variable("tr",new Range(1,1,y,x));
		
		float[][] tdata=tr.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			if(linear)
				tdata[j][i]=(float)(-1.0+2.0*j/(y-1.0)); // Linear
			else
				tdata[j][i]=(float)Math.sin(Math.PI/2.0*(-1.0+2.0*j/(y-1.0))); // Sin
		}
		
		String fname=linear?"Linear":"Sin";
		
		if(res!=1) fname+=res;
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Cartesian/"+fname+".dat");
		cdws.writeData(tr); cdws.closeFile();
	}
}
