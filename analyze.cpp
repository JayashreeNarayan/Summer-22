//Code for 2 matrices whose action is given as eqn 1 in the multi-matrx models paper and is a 0-d model, so no lattice sites, just a single point.
#include <iostream>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <random>
#include <chrono>
#include <complex>

using namespace std;

double rando()     //function to generate random numbers between [0,1) an analouge od 'drand48' but for windows.
{
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(0, 1);
    //cout << unif(rng) << endl;
    // ready to generate random numbers
    return unif(rng);
}
int main()
{   typedef complex<double> comp;
    static int first_time = 1;
    static ofstream f_corr_X, f_corr_Y, f_site_X, f_site_Y, f_a_rate_X, f_a_rate_Y, ftrace;

    if (first_time)
    {
       f_corr_X.open("corr_matrix_X.txt");      //file to save the correlator of X and Y (below)
        if (f_corr_X.bad())
        {cout << "Failed to open correlator file\n"<< flush;}
       f_corr_Y.open("corr_matrix_Y.txt");
        if (f_corr_Y.bad())
        {cout << "Failed to open correlator file\n"<< flush;}
       f_a_rate_X.open("acceptance_matrix_X.txt");      //file to save the acceptance rate of X and Y (below)
        if (f_a_rate_X.bad())
        {cout << "Failed to open acceptance rate file\n"<< flush;}
       f_a_rate_Y.open("acceptance_matrix_Y.txt");
        if (f_a_rate_Y.bad())
        {cout << "Failed to open acceptance rate file\n"<< flush;}
       f_site_X.open("site_matrix_X.txt");          //file to save the values of the site matrices of X and Y (below)
        if (f_site_X.bad())
        {cout << "Failed to open sites data file\n"<< flush;}
       f_site_Y.open("site_matrix_Y.txt");
        if (f_site_Y.bad())
        {cout << "Failed to open sites data file\n"<< flush;}
        ftrace.open("trace.txt");      //file to save the trace of all the terms in action
        if (ftrace.bad())
        {cout << "Failed to open trace file\n"<< flush;}
        first_time = 0;
    }

    int N=5;              //dimension of matrix
    int THERM = 0;        //number of thermalization steps
    int SWEEPS = 10000;   //number of generation steps
    int GAP = 100;        //interval between measurements
    double DELTA = 0.5;   //random shift range, −DELTA <= delta <= DELTA

    double u=0.0, L=0.5;         // the random number generated for the sake of probability and L is the value of Lambda
    int accept_X=0, accept_Y=0.0;     //acceptance rate of X and Y 

    double w=1.0, m=1.0, re, im;    //w is the frequency and m is the mass, dS is the change in action, re and im are varibles to get the real and imaginary parts of a complex number

    comp siteX[N][N], siteY[N][N];      //matrix at the point 
    comp shiftX[N][N], shiftY[N][N];    //randomly generated matrix which has many random numbers
    comp ositeX[N][N], ositeY[N][N];    //matrix before the shift is created
    comp nsiteX[N][N], nsiteY[N][N];    //matrix after the shift is created
    
    //double corr_sq[T], std_err[T]; // all error terms
   // double xsq = 0.0, xsq_sq = 0.0, x_val = 0.0, x_val_sq = 0.0, std_err_x_val=0.0, std_err_xsq=0.0;    //  TO BE FIGURED OUT
                                                                                                        //ALSO FIGURE OUT THE EIGENVLAUES THING 

    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
            {
                siteX[r][c]=siteY[r][c]=shiftX[r][c]=shiftY[r][c]=ositeX[r][c]=ositeY[r][c]=nsiteX[r][c]=nsiteY[r][c]=0.0;
            }

for(int s=0; s<SWEEPS; s++)
{
    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
            {
                ositeX[r][c]=siteX[r][c];
                ositeY[r][c]=siteY[r][c];
            }
    
    for (int r = 0; r < N; r++)
    {
        for (int c = r; c < N; c++)
        {
            if ((r == c) && (r != N - 1))
            {
                re = 2.0 * DELTA * (rando() - 0.5);
                im = 2.0 * DELTA * (rando() - 0.5);
                shiftX[r][c] = comp(re, im); // makes only the (0,0) element random

                re = 2.0 * DELTA * (rando() - 0.5);
                im = 2.0 * DELTA * (rando() - 0.5);
                shiftY[r][c] = comp(re, im); // makes only the (0,0) element random
                
            }
            if (r != c)
            {
                re = 2.0 * DELTA * (rando() - 0.5);
                im = 2.0 * DELTA * (rando() - 0.5);
                shiftX[r][c] = comp(re, im); // generating random numbers for the non-diagonal elements

                shiftX[c][r] = shiftX[r][c]; // To make it symmetric matrix

                re = 2.0 * DELTA * (rando() - 0.5);
                im = 2.0 * DELTA * (rando() - 0.5);
                shiftY[r][c] = comp(re, im); // generating random numbers for the non-diagonal elements

                shiftY[c][r] = shiftY[r][c]; // To make it symmetric matrix
            }
        }
    }
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
        {
            double re = real(shiftX[r][c]);
            double im = imag(shiftX[r][c]);
            if (r < c)
            {
                shiftX[r][c] = comp(re, -im); // to make the matrix hermitian wrt to the non-diagonal elements
            }
            re,im=0.0;

            re = real(shiftY[r][c]);
            im = imag(shiftY[r][c]);
            if (r < c)
            {
                shiftY[r][c] = comp(re, -im); // to make the matrix hermitian wrt to the non-diagonal elements
            }
        }
    double sum_X, sum_Y;
    for(int k=0; k<N; k++)
        {
            sum_X+=real(shiftX[k][k]);
            sum_Y+=real(shiftY[k][k]);
        }
    
    for (int k = 0; k < N; k++)
    {
        shiftX[k][k] = real(shiftX[k][k]);
        shiftY[k][k] = real(shiftY[k][k]);

    }

    sum_X=sum_X/N;
    sum_Y=sum_Y/N;

    //comp sum = (0, 0);
    for (int k = 0; k < N - 1; k++)
        {   shiftX[k][k]=real(shiftX[k][k])-sum_X;      //finally we have generated matrices that are hermitian
            shiftY[k][k]=real(shiftY[k][k])-sum_Y;
        }

    for(int r=0; r<N; r++)
           { 
           for(int c=0; c<N;c++)
                {
                    nsiteX[r][c]=siteX[r][c]+shiftX[r][c];
                    nsiteY[r][c]=siteY[r][c]+shiftY[r][c];
                }
            }

    comp X2n[N][N], Y2n[N][N], XYn[N][N], YXn[N][N], XYcommn[N][N], XY2n[N][N];
    comp X2o[N][N], Y2o[N][N], XYo[N][N], YXo[N][N], XYcommo[N][N], XY2o[N][N]; 

    double Tr1=0.0;
    double Tr2=0.0;
    double Tr3=0.0;

    //Calculating the value of S_new
    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
        {   X2n[r][c]=0.0;    //the matrices that will hold the results of all the necessary matrices
            Y2n[r][c]=0.0;
            XYn[r][c]=0.0;
            YXn[r][c]=0.0;
            
            for(int k=0; k<N; k++)
                {
                    X2n[r][c]+=nsiteX[r][k]*nsiteX[k][c];        //does X^2
                    Y2n[r][c]+=nsiteY[r][k]*nsiteY[k][c];        //does Y^2
    
                    XYn[r][c]+=nsiteX[r][k]*nsiteY[k][c];        //does XY
                    YXn[r][c]+=nsiteY[r][k]*nsiteX[k][c];        //does YX
                }

        }
    
    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
            {
                XYcommn[r][c]=XYn[r][c]-YXn[r][c];
            }

    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
        {   XY2n[r][c]=0.0;
            for(int k=0; k<N; k++)
                {
                    XY2n[r][c]=XYcommn[r][k]*XYcommn[k][c];
                }
        }

//Calculating the value of S_old
    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
        {   X2o[r][c]=0.0;    //the matrices that will hold the results of all the necessary matrices
            Y2o[r][c]=0.0;
            XYo[r][c]=0.0;
            YXo[r][c]=0.0;
            
            for(int k=0; k<N; k++)
                {
                    X2o[r][c]+=ositeX[r][k]*ositeX[k][c];        //does X^2
                    Y2o[r][c]+=ositeY[r][k]*ositeY[k][c];        //does Y^2
    
                    XYo[r][c]+=ositeX[r][k]*ositeY[k][c];        //does XY
                    YXo[r][c]+=ositeY[r][k]*ositeX[k][c];        //does YX
                }

        }
    
    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
            {
                XYcommo[r][c]=XYo[r][c]-YXo[r][c];
            }

    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
        {   XY2o[r][c]=0.0;
            for(int k=0; k<N; k++)
                {
                    XY2o[r][c]=XYcommo[r][k]*XYcommo[k][c];
                }
        }

    //calculating the trace of each term after taking dS
    for(int k=0; k<N; k++)
    {
        Tr1+=real(X2n[k][k]-X2o[k][k]);
        Tr2+=real(Y2n[k][k]-Y2o[k][k]);
        Tr3+=real(XY2n[k][k]-XY2o[k][k]);
    }

    ftrace<<Tr1<<" "<<Tr2<<" "<<Tr3<<endl;        //at every sweep we are finding the value of the trace of each of the terms in dS
    double dS=Tr1+Tr2-(L/N)*Tr3;                   //calculating the final value of dS

    double e=exp(-dS);

    u=rando();
    if(u<e)
    {
        for(int r=0; r<N; r++)
            for(int c=0; c<N; c++)
               {
                    siteX[r][c]=nsiteX[r][c];
                    siteY[r][c]=nsiteY[r][c];
               }
    accept_X++;     //SOMETHING WRONG HERE
    accept_Y++;
    }
    else
    {
        for(int r=0; r<N; r++)
            for(int c=0; c<N; c++)
            {
                siteX[r][c]=ositeX[r][c];
                siteY[r][c]=ositeY[r][c];
            }
    }

    if(s%GAP==0)
    {
        for(int r=0; r<N; r++)
            for(int c=0; c<N; c++)
                {f_site_X<<siteX[r][c]<<" ";
                    if((r==N)&&(c=N))
                    {f_site_X<<endl;}
                }

        for(int r=0; r<N; r++)
            for(int c=0; c<N; c++)
                {f_site_Y<<siteY[r][c]<<" ";
                    if((r==N)&&(c=N))
                    {f_site_Y<<endl;}
                }

    f_a_rate_X<<(double)accept_X*100/(s+1) <<endl; 
    f_a_rate_Y<<(double)accept_Y*100/(s+1) <<endl; 
    }

}//end of sweeps

}//end of main