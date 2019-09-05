//// Differential Evolution for solving multimodal functions
public class DEmain {
	public static void main (String arg[]) {
		System.out.println("DE: Differential Evolution");
		
		///////////////////////////
		// Initial configuration //
		///////////////////////////
		
		// #Benchmark number #The global minimum 
		// 1:Rosenbrock function constrained with a cubic and a line, f(1.0,1.0)=0
		// 2:Rosenbrock function constrained to a disk, f(1.0,1.0)=0
		// 3:Mishra's Bird function - constrained, f(-3.1302468,-1.5821422)=-106.7645367
		// 4:Simionescu function, f(+-0.84852813,-+0.84852813)=-0.072
		// 5:Rastrigin function, f(0, ..., 0)=0
		// 6:Rosenbrock function, f(1, ..., 1)=0
		// 7:Styblinski-Tang function, f(-2.903534, ..., -2.903534)=-39.16599*dim
		// 8:Ackley function, f(0, ..., 0)=0
		// 9:Levy function, f(1, ..., 1)=0
			
		int fnum=5;			// the function number to solve
		int dim=3;			// the number of dimension of the decision variable
		int trial=1;		// the number of trials with different random initial
		int ite=1000;		// the number of iterations for a trial
		int pop=100;		// the population size
		int rseed=1;		// random seed: MT method			
		int base=1;			// the base vector selection method, 1:rand, 2:best
		int cross=2;		// the crossover method, 1: binary, 2: exponential
		double dw=0.7;		// the differential weight [0,2]
		double cr=0.5;		// the crossover rate [0,1]
		double cp=0.90;		// the crossover probability [0,1]
		double cexp=1.002;	// exponential reduction
		
		
		
		System.out.printf("General Parameters:\n");
		System.out.printf("trial\titeration\tparticles\n");
		System.out.printf("%d\t%d\t\t%d\n",trial,ite,pop);
		System.out.printf("dw\tcr\tcp\n");
		System.out.printf("%.2f\t%.2f\t%.2f\n",dw,cr,cp);
		
		System.out.printf("Function:\n");
		ObjFunc.FuncName(fnum);
		
		// Answer check
//		double[] XX = {1, 1};
//		double ans=ObjFunc.EvalFunc(XX, fnum, dim);
//		System.out.printf("%f\n",ans);
		
		/////////////////////////////////////////////////
		// Individual Generation & Initial Evaluation  //
		/////////////////////////////////////////////////
		double[][] Frec = new double [ite][trial];			// the best fitness log
		double[][][] X1rec = new double [ite][dim][trial];	// the best variables log
		
		Sfmt rnd = new Sfmt(rseed);

		for(int l=0;l<trial;l++){
			
			double[][] X = new double[pop][dim];	// the decision variable matrix
			double[] F = new double[pop];		// the objective function for the personal best
			double[] Gb = new double[dim];		// the global best vector
			double[][] xylim=ObjFunc.GetLimit(fnum,dim);
			
			double minfnc=1000;
			int minind=0;
			
			for(int i=0;i<pop;i++) {
				for(int j=0;j<dim;j++) {
					X[i][j]=rnd.NextUnif()*(Math.abs(xylim[0][j])+Math.abs(xylim[1][j]))-Math.abs(xylim[0][j]);
				}
				F[i]=ObjFunc.EvalFunc(X[i], fnum, dim);	// variable evaluation
//				System.out.println(F[i]);
				if(minfnc>F[i]) {
					minfnc=F[i];
					minind=i;
					for(int j=0; j<dim; j++) {
						Gb[j]=X[minind][j];
					}
				}
			}
			
			
			///////////////////////////
			//  DE solution search   //
			///////////////////////////

//			System.out.println(Arrays.toString(anum));
//			System.out.println(F[anum[1]]);	// refer by this manner
			
			DEOperator DE = new DEOperator();	// DE operator class
			double cr1 = cr;
			
			for(int k=0; k<ite; k++) {		// MVMO loop start
				///// Differential Mutation /////
				double[][] Xm = new double [pop][dim];	// mutation vector
				Xm = DE.Mutation(dw, X, minind, rseed, base);
				double[][] Xc = new double [(int)(pop*cp)][dim];	// mutation vector
				
				///// Differential Crossover /////
				if(cross==2) {
					cr1 = cr / (Math.pow(cexp, k));
				}
				Xc = DE.Crossover(X, Xm, cr1, cp, rseed);
				
				
		
				/*
				System.out.printf("Iteration number:\t");
				System.out.printf("%d-",k+1);
				System.out.println(l+1);
				System.out.printf("The best fitness:\t");
				System.out.println(F1[anum[0]]);
				*/
			}
		}		
		        
		
		// write out the objective function value & search log
		String pathname = "C:\\result\\MVMO\\";
		String sufi= ".csv";
		String fnameF = "Fitness";
		
		WriteResult.Output(Frec, ite, trial, pathname + fnameF + sufi);
		
		double[][] Xrec = new double [ite][dim];
		for(int m=0; m<trial; m++){
			for(int i=0; i<ite; i++){
				for(int j=0; j<dim; j++){
					Xrec[i][j]=X1rec[i][j][m];
				}
			}
			String fnameX1 = "x" + (m+1);
			WriteResult.Output(Xrec, ite, dim, pathname + fnameX1 + sufi);
		}
		
		System.out.println(pathname + "test" + sufi); 
		
	}
}
