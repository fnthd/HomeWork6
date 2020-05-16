
public class HomeWork6 {

	
	
	/*
	 * Problem one part 1) Using a=a(x)=x-1/2 and only an initial condition u(x,0)=x(1-x)
	 * Parameters: delta t, delta x and time step 
	 * Output: will print out the (xj,tn) values for nUj depending on given timestep 
	 */
	 
	static public void UpwindScheme1(double deltaT, double deltaX, int timeStep) {
		
		double x = 0;
		double C1 = 0;
		
		
		//Finding the number of meshpoints 
		double meshPoints = 1.00/deltaX;
		long converter = Math.round(meshPoints);
		int arraySize = Math.toIntExact(converter);
		
		
		//finding the characteristic equations 
		//solution of dx/dt=a=x-1/2  ===> x=C1e^t+1/2
		
		//Array to store constant values for the meshpoints values for the charaterstic equations 
		double[] Constant_aj = new double [arraySize+1];
		
		for(int i=0; i<Constant_aj.length; i++) {
			C1 = x -.5; 
			Constant_aj[i] = C1;
			x = x+deltaX;
			//System.out.println(Constant_aj[i]);
		}
		
		//Setting up inital U^0_j array  
		double[] U0 = new double [arraySize+1];
		
		
		
		//X(0)  
		double X_0 = 0;
		//U(0)  
		double U_0 = 0;
		//System.out.println();
		//System.out.println();
		
		
		
		
		//Loop to fill find all the U_0_j values using the array of constat in aj
		for(int i=0; i<U0.length; i++) {
			//finding X(0)=xj where X(t)= C1*e^t+1/2
			X_0 = Constant_aj[i]*Math.pow(Math.E, 0)+.5;
			//finding U(0)=f(xj) U(t)=C C=f(xj) where f(x)=x(1-x)
			U_0 = X_0*(1-X_0);
			//saving U^0_j value ot the array 
			U0[i] = U_0;
			
			//System.out.println(U0[i]);
		
		}
		
		//System.out.println();
		//System.out.println();
		
		
		//Now we want to the value of a(xj) at each xj 
		double[] Aj = new double [arraySize+1];
		x = 0; 
		
		
		for(int i=0; i<Aj.length; i++) {
		Aj[i] = x-(1.00/2.00);
		x = x+deltaX;
		//System.out.println(Aj[i]);
		}
		
		//System.out.println();
		//System.out.println();
		
		
	    ///////////////////////////////////////////////////////Now we can apply the upwind scheme//////////////////////////////////////////////////////////////////// 
		
		//Creating the an array to store U^n+1_j 
		double[] U1 = new double [arraySize+1];
		
		
		
		int counter = 1;
		
		while (counter < timeStep+1) {
		//System.out.println("TimeStep: "+ counter );
		
		//x = 0;
		
		//for loop to find the U^n+1_j values 
		for(int i=0; i<Aj.length; i++) {
			
		//checking to see if aj is postive or negative 	
		if(Aj[i]<0) {
		
		//apply upwind scheme 	
		U1[i]= U0[i]-((Aj[i]*(deltaT/deltaX))*(U0[i+1]-U0[i]));
		}
		else {
		U1[i] = U0[i]-((Aj[i]*(deltaT/deltaX))*(U0[i]-U0[i-1]));	
		}
		
		//System.out.println("(" + x + "," + U1[i] + ")");
		//x = x+deltaX;
		}
		
		
		//System.out.println();
		//System.out.println();
		
		//Setting the U1 array to U0 to do next time step.
		for(int i=0; i<U1.length; i++) {
			
			U0[i] = U1[i];
			
			
		}
		
		
		
		counter++;
		}
		
		x=0;
		System.out.println("TimeStep: "+ timeStep);
		for(int i=0; i<U0.length; i++) {
		System.out.println("(" + x + "," + U0[i] + ")");
		x =x+deltaX;
		}
		
	}
	
	
	
	
	/*
	 * Problem one part 2) using a(x)=1/2-x with boundary conditions U(0,t)=U(1,t)=0
	 * Parameters: delta t, delta x and time step 
	 * Output: will print out the (xj,tn) values for nUj depending on given timestep
	 */
static public void UpwindScheme2(double deltaT, double deltaX, int timeStep) {
		
		double x = 0;
		double C1 = 0;
		
		
		//Finding the number of meshpoints 
		double meshPoints = 1.00/deltaX;
		long converter = Math.round(meshPoints);
		int arraySize = Math.toIntExact(converter);
		
		
		//finding the characteristic equations 
		//solution of dx/dt=a=1/2-x ====> x=C1e^-t+1/2
		
		//Array to store constant values for the meshpoints values for the charaterstic equations 
		double[] Constant_aj = new double [arraySize+1];
		
		
		
		for(int i=0; i<Constant_aj.length; i++) {
			C1 = x-(1.00/2.00); 
			Constant_aj[i] = C1;
			x = x+deltaX;
			
		}
		
		//Setting up intial U^0_j array  
		double[] U0 = new double [arraySize+1];
		
		
		
		//X(0)  
		double X_0 = 0;
		//U(0)  
		double U_0 = 0;
		
		
		//Applying the boundary conditions to U(0,1)=0 
		U0[0] = 0; 
		
		
		//Loop to fill find all the U_0_j values using the array of constant in aj until Uj-1 (Uj will be set to zero by boundary condition) 
		for(int i=1; i<U0.length-1; i++) {
			//finding X(0)=xj where X(t)= C1*e^t+1/2
			X_0 = Constant_aj[i]*Math.pow(Math.E, 0)+.5;
			//finding U(0)=f(xj) U(t)=C C=f(xj) where f(x)=x(1-x)
			U_0 = X_0*(1-X_0);
			//saving U^0_j value ot the array 
			U0[i] = U_0;
			
			
		
		}
		
		
		//Applying the boundary conditions to U(1,0)=0
		U0[arraySize] = 0;
		
		
		
		System.out.println();
		System.out.println();
		
		
		//Now we want to the value of a(xj) at each xj 
		double[] Aj = new double [arraySize+1];
		x = 0; 
		
		
		for(int i=0; i<Aj.length; i++) {
		Aj[i] = x-(1.00/2.00);
		x = x+deltaX;
		//System.out.println(Aj[i]);
		}
		
		//System.out.println();
		//System.out.println();
		
		
	    ///////////////////////////////////////////////////////Now we can apply the upwind scheme//////////////////////////////////////////////////////////////////// 
		
		//Creating the an array to store U^n+1_j 
		double[] U1 = new double [arraySize+1];
		
		
		
		int counter = 1;
		
		while (counter < timeStep+1) {
		//System.out.println("TimeStep: "+ counter );
		
		//x = 0;
			
			
			
		//assigning U(0,t+1) = 0	for boundary condition 
		U1[0] = 0;	
			
			
			
		
		//for loop to find the U^n+1_j values 
		for(int i=1; i<Aj.length-1; i++) {
			
		//checking to see if aj is postive or negative 	
		if(Aj[i]<0) {
		
		//apply upwind scheme 	
		U1[i]= U0[i]-((Aj[i]*(deltaT/deltaX))*(U0[i+1]-U0[i]));
		}
		else {
		U1[i] = U0[i]-((Aj[i]*(deltaT/deltaX))*(U0[i]-U0[i-1]));	
		}
		
		//System.out.println("(" + x + "," + U1[i] + ")");
		//x = x+deltaX;
		}
		
		//Setting up the last value of U(1,t+1) =0
		U1[0]=0;
		
		
		
		//System.out.println();
		//System.out.println();
		
		//Setting the U1 array to U0 to do next time step.
		for(int i=0; i<U1.length; i++) {
			
			U0[i] = U1[i];
			
			
		}
		
		
		
		counter++;
		}
		
		
		
		x=0;
		System.out.println("TimeStep: "+ timeStep);
		for(int i=0; i<U0.length; i++) {
		System.out.println("(" + x + "," + U0[i] + ")");
		x =x+deltaX;
		}
	}
	

        ////////////////////////////////////////////end of methods for problem one 


        /*
         * Problem 2 code 
         * Finite difference method to solve the poisson equation  equation 
         */




static public void FiniteDifference(double deltaY, double deltaX) {
	
	//Solving the Poisson equation uxx+uyy+5pi^2(sin(xpi)sin(2ypi)) with 0 Dirichlet boundary conditions  
	//Finding J where J= 1/deltaX 
			double meshPoints = 1.00/deltaX;
			long converter = Math.round(meshPoints);
			int J2 = Math.toIntExact(converter);
			int J = J2*2; 
			int arraySize1 = J*2+1;
			int arraySize = (J+1)*(J+1);
			
			//Q=-4/(deltaX*deltaY)
			double Q = -4.00/(deltaY*deltaX);
			//R = 1/(deltaX*deltaY)
			double R = 1/(deltaX*deltaY);
			
			
			//Array to hold the U(xj,yj) where 0<=j<=J 
			double[] U = new double [arraySize];
			
			//Applying the Dirichlet boundary conditions 
			
			
			//Setting the first U(0,yj)=0 and U(xJ,yj)=0 
			for(int i=0; i<J+1; i++) {
				U[i] =0;
				U[arraySize-i-1] = 0;
			}
			
			//We finish setting up the other boundary values to 0
			int counter =J+1;
			
			for(int i=1; i<J+1; i++) {
				
				U[counter*i] =0;
				U[counter*i-1] = 0;	
			}
			
			
			
			
			
			
			//Now we build our system of equations for U(r+1,s)+U(r-1,s)+U(r,s+1)+U(r,s-1)-4U(r,s) = f(r,s) 
			//AU=f system of equations 
			
			//Setting up a matrix containing  the coefficients 
			double[][] A = new double [(J-1)*(J-1)][(J-1)*(J-1)];
			//Setting up matrix containing U(xj,yj) values 
			double[] U1 =  new double [(J-1)*(J-1)];
			//Setting up matrix containg the F(xj,yj) values 
			double[] f =  new double [(J-1)*(J-1)];
			
			//filling out the the A matrix 
			//Filling out the penta diagonal system  
			
			counter = 0; //nth row  
			int tracker = 0; //ith column 
			int tracker2 = 0; //ith column  
			int counter2 =0; //nth row  
			int tracker3 =0; //ith column
			int counter3 =0; //nth row 
			
			
		//For loop to fill 	
		for(int j = 1; j<J; j++) {
				
			for(int i = 1; i<J-1; i++) {
				
				
				A[counter][tracker] =Q;  //Main diagonal
				tracker++; //next column 
				A[counter][tracker] = R; //upper diagonal 
				counter++; //next row 
				A[counter][tracker-1]= R; //lower diagonal 
				
			
			}
			   
			
			   A[counter][tracker] =Q;
			  
			   
			   counter++;
			   tracker++;
			   //setting up the lower matrix of A
			   if(j<J-1) {
				   counter2=counter;
				   for(int k=0; k<J-1; k++) {
				   A[counter2][tracker2] =R;
				   tracker2++;
				   counter2++;
				   
				   }
			   }
			   
			   //setting up the upper matrix of A
			   if(j<J-1) {
				   tracker3=tracker;
				   for(int k=0; k<J-1; k++) {
				   A[counter3][tracker3] =R;
				   tracker3++;
				   counter3++;
				   
				   }
			   }   
		}
			 
			
			
			System.out.println("Matrix A to be solved");
			
			for(int i =0; i<A.length;i++) {
				
				for(int j =0; j<A.length; j++) {
					System.out.print(A[i][j]+ " ");
					
				}
				System.out.println();
				
			}
			
			double x =-1+deltaX;
			double y =-1+deltaY;
			int index =0;
			//Setting up the values of f(xj,yj) 
			for(int i=0; i<J-1; i++) {
				
				for(int j =0; j<J-1; j++) {
					f[index] = -5*Math.pow(Math.PI, 2)*Math.sin(x*Math.PI)*Math.sin(2*y*Math.PI);
					y = y+deltaY;
					index++;
					
				}
				 x=x+deltaX; 
				 y=-1+deltaY;
			}
			
			
			
			//Matrix A is setup and f now we solve for U, AU=f by calling Matrix solver where guassian elimination will be used to find the U values  
			double[] U0 = MatrixSolver(A,f);
			
			
			
			
			
			
			int index2 =J+2;
			int index3 = 1;
			
			for(int i =1; i<J-1;i++){
				
				for(int j=1; j<J;j++) {
				U[index2]=U0[index3];
				index2++;
				index3++;
				
				}
				index2++;
				index2++;
				
			}
			
			//Setting up the final array 
			for(int i=0;i<J-2;i++) {
				U[index2]=U0[index3]; 
				index2++;
				index3++;
			}
			
			 x =-1;
			 y =-1;
			 index =0;
			 

			System.out.println("Final solution");
			for(int i=0; i<J+1;i++) {
				for(int j = 0; j<J+1; j++) {
					System.out.println("("+ x + ", "+y+", "+ U[index]+ ")");
					y = y+deltaY;
					index++;
				}
				x=x+deltaX; 
				y=-1;
 			}
}



//Matrix solver for the finite difference method 
static public double checker = 1e-10;
public static double[] MatrixSolver(double[][] A, double[] b) {
    int n = b.length;

    for (int p = 0; p < n; p++) {

        // find pivot row and swap
        int max = p;
        for (int i = p + 1; i < n; i++) {
            if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                max = i;
            }
        }
        double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
        double   t    = b[p]; b[p] = b[max]; b[max] = t;

        // singular or nearly singular
        if (Math.abs(A[p][p]) <= checker) {
            throw new ArithmeticException("Matrix is singular");
        }

        // pivot within A and b
        for (int i = p + 1; i < n; i++) {
            double alpha = A[i][p] / A[p][p];
            b[i] -= alpha * b[p];
            for (int j = p; j < n; j++) {
                A[i][j] -= alpha * A[p][j];
            }
        }
    }

    // back substitution
    double[] x = new double[n];
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
   
    
    return x;
}






static public void FiniteDifference2(double deltaY, double deltaX) {
	
	
	double meshPoints = 1.00/deltaX;
	long converter = Math.round(meshPoints);
	int J2 = Math.toIntExact(converter);
	int J = J2*2; 
	int arraySize1 = J*2+1;
	int arraySize = (J+1)*(J+1);
	
	//Q=-4/(deltaX*deltaY)
	double Q = -4.00/(deltaY*deltaX);
	//R = 1/(deltaX*deltaY)
	double R = 1/(deltaX*deltaY);
	
	//Array to hold the U(xj,yj) where 0<=j<=J 
	double[] U = new double [arraySize];
	
	
	
	//Applying boundary conditions 
	
	
	//Setting the first U(0,yj)=0 and U(xJ,yj)=2 
	for(int i=0; i<J+1; i++) {
		U[i] =0;
		U[arraySize-i-1] = 2;
	}
	
	//We finish setting up the other boundary values U(xj,0)=U(xj,J)=x+1
	int counter =J+1;
	
	double deltaX0 = -1+deltaX;
	
	for(int i=1; i<J+1; i++) {
		if(i<2) {
			U[counter*i-1] = 0;
			U[counter*i] =deltaX0+1;
			
		}
		else if(i>J-1) {
			U[counter*i-1] = deltaX0+1;
			U[counter*i] =2 ;
		}
		else {
			U[counter*i-1] = deltaX0+1;
			deltaX0 = deltaX0 +deltaX;	
			U[counter*i] =deltaX0+1;
		}
		
		
	}
	
	
	
	//Now we build our system of equations for U(r+1,s)+U(r-1,s)+U(r,s+1)+U(r,s-1)-4U(r,s) = 0 
	//AU=f system of equations 
	
	//Setting up a matrix containing  the coefficients 
	double[][] A = new double [(J-1)*(J-1)][(J-1)*(J-1)];
	//Setting up matrix containing U(xj,yj) values 
	double[] U1 =  new double [(J-1)*(J-1)];
	//Setting up matrix containg the F(xj,yj) values 
	double[] f =  new double [(J-1)*(J-1)];
	
	//filling out the the A matrix 
	//Filling out the penta diagonal system  
	
	counter = 0; //nth row  
	int tracker = 0; //ith column 
	int tracker2 = 0; //ith column  
	int counter2 =0; //nth row  
	int tracker3 =0; //ith column
	int counter3 =0; //nth row 
	
	
//For loop to fill 	
for(int j = 1; j<J; j++) {
		
	for(int i = 1; i<J-1; i++) {
		
		
		A[counter][tracker] =Q;  //Main diagonal
		tracker++; //next column 
		A[counter][tracker] = R; //upper diagonal 
		counter++; //next row 
		A[counter][tracker-1]= R; //lower diagonal 
		
	
	}
	   
	
	   A[counter][tracker] =Q;
	  
	   
	   counter++;
	   tracker++;
	   //setting up the lower matrix of A
	   if(j<J-1) {
		   counter2=counter;
		   for(int k=0; k<J-1; k++) {
		   A[counter2][tracker2] =R;
		   tracker2++;
		   counter2++;
		   
		   }
	   }
	   
	   //setting up the upper matrix of A
	   if(j<J-1) {
		   tracker3=tracker;
		   for(int k=0; k<J-1; k++) {
		   A[counter3][tracker3] =R;
		   tracker3++;
		   counter3++;
		   
		   }
	   }   
}
	 
	
	
	System.out.println("Matrix A to be solved");
	
	for(int i =0; i<A.length;i++) {
		
		for(int j =0; j<A.length; j++) {
			System.out.print(A[i][j]+ " ");
			
		}
		System.out.println();
		
	}
	
	double x =-1+deltaX;
	double y =-1+deltaY;
	int index =0;
	//Setting up the values of f(xj,yj) =0
	for(int i=0; i<J-1; i++) {
		
		for(int j =0; j<J-1; j++) {
			f[index] = 0;
			y = y+deltaY;
			index++;
			
		}
		 x=x+deltaX; 
		 y=-1+deltaY;
	}
	
	
	//Matrix A is setup and f now we solve for U, AU=f by calling Matrix solver where guassian elimination will be used to find the U values  
	double[] U0 = MatrixSolver(A,f);
	
	
	int index2 =J+2;
	int index3 = 1;
	
	for(int i =1; i<J-1;i++){
		
		for(int j=1; j<J;j++) {
		U[index2]=U0[index3];
		index2++;
		index3++;
		
		}
		index2++;
		index2++;
		
	}
	
	
	 x =-1;
	 y =-1;
	 index =0;
	
	
	System.out.println("Final solution");
	for(int i=0; i<J+1;i++) {
		for(int j = 0; j<J+1; j++) {
			System.out.println("("+ x + ", "+y+", "+ U[index]+ ")");
			y = y+deltaY;
			index++;
		}
		x=x+deltaX; 
		y=-1;
		}
	
	
	
}



	
//////////////////////////////////////////////////////////////////////////////////Main
	
	public static void main(String[] args) {
		
		
		//info for problem one 
		double deltax =.1;
		double deltat =.001;
		int timestep = 15;
		
		
		System.out.println("Applying the upwind scheme for problem 1 part 1"); 
		System.out.println("a=x-1/2");
		System.out.println("Initial Condition U(X,0)=x(1-x)");
		System.out.println("Print out U(xj,tn)");
		System.out.println("(xj,tn)");
		System.out.println();
		System.out.println();
		UpwindScheme1(deltat,deltax,timestep);
		
		System.out.println();
		System.out.println();
		
		System.out.println("Applying the upwind scheme for problem 1 part 2"); 
		System.out.println("a=1/2-x");
		System.out.println("Initial Condition U(X,0)=x(1-x)");
		System.out.println("Boundary conditions: U(0,t)=U(1,t)=0");
		System.out.println("Print out U(xj,tn)");
		System.out.println("(xj,tn)");
		UpwindScheme2(.001,.1,15);
		
		System.out.println();
		System.out.println();
		
		//info for finite difference 
		double deltax2 =1.00/3.00;
		double deltay = deltax2;
		
		
		System.out.println("Applying the finite difference method for problem 2"); 
		System.out.println("To solve Poisson Equation uxx+uyy+5pi^2(sin(xpi)sin(2ypi)");
		System.out.println("With 0 Dirichlet boundary conditon");
		System.out.println("solution will print out the points (xj,yj,zj)");
		FiniteDifference(deltay,deltax2);
		
		System.out.println();
		System.out.println();
		
		System.out.println("Applying the finite difference method for problem 3 ");
		System.out.println("Boundary conditions U(x,y)==0 if x=-1, ==2 if x=1, ==x+1 if y=-1 or y=1");
		System.out.println("solution will print out the points (xj,yj,zj)");
		FiniteDifference2(deltay,deltax2);
        
		
	}

}
