//
import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;

import diffuse.DiffusionModel;
import miniufo.application.statisticsModel.LagrangianStatisticsByDavis;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.lagrangian.LSM2nd;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.lagrangian.StochasticModel;
import miniufo.lagrangian.StochasticParams;
import miniufo.lagrangian.StochasticModel.BCType;
import miniufo.util.Region2D;

//
public final class StochasticModels{
	// tracking parameters
	static final int  intLen=730;		// length of integration
	static final float DiffX=10000;		// m^2/s
	static final float DiffY=10000;		// m^2/s
	static final float TvX  =4;			// days
	static final float TvY  =4;			// days
	static final float TaX  =4;			// days
	static final float TaY  =4;			// days
	
	static final float[][] diff=new float[][]{{DiffX,0},{0,DiffY}};
	static final float[][] Tvel=new float[][]{{TvX,1},{1,TvY}};
	static final float[][] Tacc=new float[][]{{TaX,1},{1,TaY}};
	
	// general parameters
	static final boolean writeTraj=false;
	
	static final String EulerCTL=
		"dset ^ctl.dat\n"+
		"title template\n"+
		"undef -9999\n"+
		"xdef 31 linear 140 2\n"+
		"ydef 11 linear -10 2\n"+
		"zdef  1 linear   1 1\n"+
		"tdef  1 linear 1Jan2001 1dy\n"+
		"vars 1\n"+
		"u 0 99 u\n"+
		"endvars\n";
	
	static final String path="/lustre/home/qianyk/Data/Idealized/Validate/";
	
	
	/** test*/
	public static void main(String[] args){
		generateLagrangianData("StaticMean");
		cLagrangianStatistics(new Region2D(168f,-2,172,2),180,"StaticMean");
	}
	
	/**
	 * generate synthetic Lagrangian drifter data
	 * 
	 * @param ctlname	mean flow data
	 * @param ensemble	ensemble number
	 * @param intLen	length of integration (steps)
	 * 
	 * @return	ps		simulated particles
	 */
	static void generateLagrangianData(String tag){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+tag+".ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Region2D wholeregion=dd.toRegion2D();
		
		Function<Record,StochasticParams> f1=r->{
			//float lon=r.getLon();
			
			//float Diff=DiffX+(10-1)*DiffX*(lon-140)/(160-140);
			return new StochasticParams(Tvel,Tacc,diff);
		};
		
		StochasticModel sm=new LSM2nd(1,dd,BCType.Reflected,BCType.Reflected,f1);
		sm.setVelocityBuffer("u","v",1);	// set initial velocity buffer
		
		List<Particle> ps=sm.deployPatch(wholeregion,0.5f,1,intLen);
		
		System.out.println("\nrelease "+ps.size()+" particles\n");
		
		sm.simulateParticles(ps,"u","v",intLen);
		
		System.out.println("finished");
		
		if(writeTraj) DiffusionModel.writeTrajAndGS(ps,path+"TXT/",wholeregion);
		
		DiffusionModel.writeParticleList(path+tag+"LD.dat",ps);
	}
	
	static void cLagrangianStatistics(Region2D region,int tRad,String tag){
		List<Particle> ps=DiffusionModel.readParticleList(path+tag+"LD.dat");
		
		DataDescriptor dd=DiagnosisFactory.parseContent(EulerCTL).getDataDescriptor();
		
		//EulerianStatistics estat=new EulerianStatistics(ps,dd,false);
		//estat.removeMeansOfBins();
		
		/**************** Lagrangian statistics ****************/
		System.out.println("\nLagrangian Statistics...");
		LagrangianStatisticsByDavis lstat=new LagrangianStatisticsByDavis(ps,dd);
		
		Predicate<Record> cond=r->region.inRange(r.getXPos(),r.getYPos());
		
		lstat.cStatisticsByDavisTheory     (cond,tRad).toFile(path+"Diff/Lstat"+tag+"1.txt");
		lstat.cStatisticsByTaylorTheory    (cond,tRad).toFile(path+"Diff/Lstat"+tag+"2.txt");
		lstat.cStatisticsByDispersionTheory(cond,tRad).toFile(path+"Diff/Lstat"+tag+"3.txt");
	}
}
