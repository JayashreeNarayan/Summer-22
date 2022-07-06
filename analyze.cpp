//Code for 2 matrices whose action is given as eqn 1 in the multi-matrx models paper and is a 0-d model, so no lattice sites, just a single point.
#include <iostream>
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
    static ofstream f_a_rate, f_trace;

    if (first_time)
    {
       f_a_rate.open("data/acceptance_rate.txt");      //file to save the acceptance rate of X and Y (below)
        if (f_a_rate.bad())
        {cout << "Failed to open acceptance rate file\n"<< flush;}
        f_trace.open("data/trace.txt");      //file to save the trace of all the terms in action
        if (f_trace.bad())
        {cout << "Failed to open trace file\n"<< flush;}
        first_time = 0;
    }

    int N=5;              //dimension of matrix
    int THERM = 0;        //number of thermalization steps
    int SWEEPS = 1000000;   //number of generation steps
    int GAP = 1000;        //interval between measurements
    double DELTA = 0.07;   //random shift range, −DELTA <= delta <= DELTA

    double u=0.0, L=8.0;         // the random number generated for the sake of probability and L is the value of Lambda
    int accept_count=0 ;     //acceptance rate of X and Y 

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
        for (int c = 0; c < N; c++)
        {
            if (r == c)
            {
                re = 2.0 * DELTA * (rando() - 0.5);
                im = 0.0;
                shiftX[r][c] = comp(re, im); // makes only the (n,n) element random, diagonal elemets only real

                re = 2.0 * DELTA * (rando() - 0.5);
                im = 0.0;
                shiftY[r][c] = comp(re, im); // makes only the (n,n) element random
                
            }
            if (r > c)
            {
                re = 2.0 * DELTA * (rando() - 0.5);
                im = 2.0 * DELTA * (rando() - 0.5);
                shiftX[r][c] = comp(re, im); // generating random numbers for the non-diagonal elements

                shiftX[c][r] = comp(re, (-1.0)*im); // To make it hermitian matrix

                re = 2.0 * DELTA * (rando() - 0.5);
                im = 2.0 * DELTA * (rando() - 0.5);
                shiftY[r][c] = comp(re, im); // generating random numbers for the non-diagonal elements

                shiftY[c][r] = comp(re, (-1.0)*im); // To make it hermitian matrix
            }
        }
    }
    
    double sum_X, sum_Y;
    for(int k=0; k<N; k++)
        {
            sum_X+=real(shiftX[k][k]);
            sum_Y+=real(shiftY[k][k]);
        }

    sum_X=sum_X/N;
    sum_Y=sum_Y/N;

    //comp sum = (0, 0);
    for (int k = 0; k < N ; k++)
        {   shiftX[k][k]=comp(real(shiftX[k][k])-sum_X,0);      //finally we have generated matrices that are hermitian
            shiftY[k][k]=comp(real(shiftY[k][k])-sum_Y,0);
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
    comp X2[N][N], Y2[N][N], XY[N][N], YX[N][N], XYcomm[N][N], XYcomm2[N][N],X2Y2[N][N],XYXY[N][N]; 

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
                    XY2n[r][c]+=XYcommn[r][k]*XYcommn[k][c];
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
                    XY2o[r][c]+=XYcommo[r][k]*XYcommo[k][c];
                }
        }

    //calculating the trace of each term after taking dS
    for(int k=0; k<N; k++)
    {
        Tr1+=real(X2n[k][k]-X2o[k][k]);
        Tr2+=real(Y2n[k][k]-Y2o[k][k]);
        Tr3+=real(XY2n[k][k]-XY2o[k][k]);
    }

    //ftrace<<Tr1<<" "<<Tr2<<" "<<Tr3<<endl;        //at every sweep we are finding the value of the trace of each of the terms in dS
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
    accept_count++;   
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
    double Tr1=0.0;
    double Tr2=0.0;
    double Tr3=0.0;
    double Tr4=0.0;
    double Tr5=0.0;

    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
        {   X2[r][c]=0.0;    //the matrices that will hold the results of all the necessary matrices
            Y2[r][c]=0.0;
            XY[r][c]=0.0;
            YX[r][c]=0.0;
            
            for(int k=0; k<N; k++)
                {
                    X2[r][c]+=siteX[r][k]*siteX[k][c];        //does X^2
                    Y2[r][c]+=siteY[r][k]*siteY[k][c];        //does Y^2
    
                    XY[r][c]+=siteX[r][k]*siteY[k][c];        //does XY
                    YX[r][c]+=siteY[r][k]*siteX[k][c];        //does YX
                }

        }
    
    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
            {
                XYcomm[r][c]=XY[r][c]-YX[r][c];
                
            }

    for(int r=0; r<N; r++)
        for(int c=0; c<N; c++)
        {   XYcomm2[r][c]=0.0;
            XYXY[r][c]=0.0;
            X2Y2[r][c]=0.0;
            for(int k=0; k<N; k++)
                { 
                    XYcomm2[r][c]+=XYcomm[r][k]*XYcomm[k][c];
                    XYXY[r][c]+=XY[r][k]*XY[k][c];
                    X2Y2[r][c]+=X2[r][k]*Y2[k][c];
                }
        }


    for(int k=0; k<N; k++)
    {
        Tr1+=real(X2[k][k]);
        Tr2+=real(Y2[k][k]);
        Tr3+=real(XYcomm2[k][k]);
        Tr4+=real(XYXY[k][k]);
        Tr5+=real(X2Y2[k][k]);
    }

    f_trace<<Tr1<<" "<<Tr2<<" "<<Tr3<<" "<<Tr4<<" "<<Tr5<<endl;        //at every sweep we are finding the value of the trace of each of the terms in dS

    f_a_rate<<(double)accept_count*100/(s+1) <<endl;  
    }

}//end of sweeps

}//end of main
