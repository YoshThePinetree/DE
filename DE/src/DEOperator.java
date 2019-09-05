// DE Operator Class
public class DEOperator {
	
	// Method Mutation
	public double[][] Mutation(double F, double X[][], int rseed){

		Sfmt rnd = new Sfmt(rseed);		
		
		int dim = X[0].length;	
		int pop = X.length;
		double[][] Xm = new double[pop][dim];
		
		int r1, r2, r3;
		
		for(int i=0; i<pop; i++) {
			r1=rnd.NextInt(pop-1);	// the first mutation vector
			while(r1==i) {
				r1=rnd.NextInt(pop-1);
			}
			r2=rnd.NextInt(pop-1);	// the second mutation vector
			while(r2==i || r2==r1) {
					r2=rnd.NextInt(pop-1);
			}
			r3=rnd.NextInt(pop-1);	// the second mutation vector
			while(r3==i || r3==r1 || r3==r2) {
					r3=rnd.NextInt(pop-1);
			}
			
			for(int j=0; j<dim; j++) {
				Xm[i][j] = X[r1][j] + F*(X[r2][j]-X[r3][j]);
			}
			
		}
		
		return Xm;
	}

}
