
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in a set of x and y values and computes the line of best fit.                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void least_squares_regression(int n,double x[],double y[],double *a,double *b,double *r2)
{
    int i        = 0;                 //Standard variable used in loops
    int j        = 0;                 //Standard variable used in loops
    int k        = 0;                 //Standard variable used in loops
    double xsum  = 0;                 //Sum of x
    double ysum  = 0;                 //Sum of y
    double x2sum = 0;                 //Sum of x squared
    double y2sum = 0;                 //Sum of y squared
    double xysum = 0;                 //Sum of x*y
    double r     = 0;                 //Correlation coeficient

    for (i=0;i<n;i++) //loop over items in arrays
    {
        xsum  = xsum + x[i];          //calculate sum(x)
        ysum  = ysum + y[i];          //calculate sum(y)
        x2sum = x2sum + pow(x[i],2);  //calculate sum(x^2)
        y2sum = y2sum + pow(y[i],2);  //calculate sum(y^2)
        xysum = xysum + x[i]*y[i];    //calculate sum(x*y)
    }

    double denominator = n*x2sum - xsum*xsum;

    *a  = (n*xysum - xsum*ysum)/denominator;                                                  //calculate slope
    *b  = (x2sum*ysum - xsum*xysum)/denominator;                                              //calculate intercept
    r   = (xysum - (xsum*ysum/n)) /  sqrt((x2sum - pow(xsum,2)/n) * (y2sum - pow(ysum,2)/n)); //calculate correlation coeficient
    *r2 = r*r;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Do a nth order polynomial fit                                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d poly_fit(dv1d &data_y,dv1d &data_x,int order)
{
    int i,j,k,n,N;

    int coef_size = order+1;
    dv1d coef(coef_size,0.0);

    //set the number of data points to fit to
    N = data_x.size();

    //create array to hold data
    double x[N];
    double y[N];
    for(i=0; i<N; i++)
    {
        x[i] = data_x[i];
        y[i] = data_y[i];
       //printf("x %f y %f \n",data_x[i],data_y[i]);
    }

    //set the degree of the polynomial (how many terms)
    n = coef_size-1;

    double X[2*n+1];                            //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    for(i=0; i<2*n+1; i++)
    {
        X[i] = 0;
        for(j=0;j<N;j++)
        {
            X[i] = X[i] + pow(x[j],i);          //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
        }
    }

    double B[n+1][n+2];                         //B is the Normal matrix(augmented) that will store the equations,
    double a[n+1];                              //a is for value of the final coefficients
    for(i=0; i<=n; i++)
    {
        for(j=0; j<=n; j++)
        {
            B[i][j] = X[i+j];                   //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
        }
    }

    double Y[n+1];                              //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    for(i=0; i<n+1; i++)
    {
        Y[i] = 0;
        for(j=0; j<N; j++)
        {
            Y[i] = Y[i] + pow(x[j],i)*y[j];     //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
        }
    }

    for (i=0; i<=n; i++)
    {
        B[i][n+1] = Y[i];                       //load the values of Y as the last column of B(Normal Matrix but augmented)
    }
    n=n+1;                                      //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations

    for(i=0; i<n; i++)                          //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
    {
        for(k=i+1; k<n; k++)
        {
            if(B[i][i] < B[k][i])
            {
                for(j=0; j<=n; j++)
                {
                    double temp=B[i][j];
                    B[i][j] = B[k][j];
                    B[k][j] = temp;
                }
            }
        }
    }

    for(i=0; i<n-1; i++)                        //loop to perform the gauss elimination
    {
        for(k=i+1; k<n; k++)
        {
            double t = B[k][i]/B[i][i];
            for (j=0; j<=n; j++)
            {
                B[k][j] = B[k][j] - t*B[i][j]; //make the elements below the pivot elements equal to zero or elimnate the variables
            }
        }
    }

    for(i=n-1; i>=0; i--)                      //back-substitution
    {                                          //x is an array whose values correspond to the values of x,y,z..
        a[i] = B[i][n];                        //make the variable to be calculated equal to the rhs of the last equation
        for(j=0; j<n; j++)
        {
            if(j != i)                         //then subtract all the lhs values except the coefficient of the variable whose value is being calculated
            {
                a[i] = a[i] - B[i][j]*a[j];
            }
        }
        a[i] = a[i]/B[i][i];                   //now finally divide the rhs by the coefficient of the variable to be calculated
    }

    for(i=0; i<n; i++)
    {
        coef[i] = a[i];
        //printf("coef[%d] %f \n",i,coef[i]);
    }

    return coef;
}

