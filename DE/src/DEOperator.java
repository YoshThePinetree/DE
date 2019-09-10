// DE Operator Class
public class DEOperator {
	
	// Method Mutation (randor best)
	public double[][] Mutation(double F, double X[][], int mi, int rseed, int base){

		Sfmt rnd = new Sfmt(rseed);		
		
		int dim = X[0].length;	
		int pop = X.length;
		double[][] Xm = new double[pop][dim];
		
		int r1, r2, r3;
		
		for(int i=0; i<pop; i++) {
			
			r1=rnd.NextInt(pop-1);	// the base mutation vector for rand
			while(r1==i) {
				r1=rnd.NextInt(pop-1);
			}
			
			if(base==2) {	// if base=best: then r1 is mi
				r1=mi;
			}
			
			r2=rnd.NextInt(pop-1);	// the first differential mutation vector
			while(r2==i || r2==r1) {
					r2=rnd.NextInt(pop-1);
			}
			r3=rnd.NextInt(pop-1);	// the second differential mutation vector
			while(r3==i || r3==r1 || r3==r2) {
					r3=rnd.NextInt(pop-1);
			}
			
			for(int j=0; j<dim; j++) {
				Xm[i][j] = X[r1][j] + F*(X[r2][j]-X[r3][j]);
			}

		}
		return Xm;
	}
	
	public double[][] Crossover(double X[][], double V[][], double cr, double cp, int rseed){
		
		Sfmt rnd = new Sfmt(rseed);
		
		int dim = X[0].length;	
		int pop = X.length;
		double[][] U = new double[(int)(pop*cp)][dim];
		int[] indzero = new int[pop];	// 0-1 index for crossover
		
		// select for individuals which is not crossovered
		int stop=1;
		int m, sum=0;
		while(stop==1) {
			m = rnd.NextInt(pop);
			indzero[m]=1;
			sum = DEOperator.isum(indzero);
			if(sum >= (pop-(pop*cp))) {
				break;
			}
		}
		
		
		// Crossover
		int n=0;	// counter for U
		int jrand=0;	// randomly selected dimension
		double r;
		for(int i=0; i<pop; i++) {
			if(indzero[i]==0) {		// crossover only indzero is 0
				jrand=rnd.NextInt(dim);
				
				for(int j=0; j<dim; j++) {
					r=rnd.NextUnif();
					if(j!=jrand) {
						if(r<=cr) {
							U[n][j]=V[i][j];
						}else if(r>cr){
							U[n][j]=X[i][j];
						}						
					}else{
						U[n][j]=V[i][j];
					}
				}
				
				n++;
			}			
		}

		return U;
	}
	
	public static int isum(int A[]) {
		int sum=0;
		int len=A.length;
		for(int i=0; i<len; i++) {
			sum = sum + A[i];
		}
		return sum;
	}

	
}
