/*
    sampleHRnew.cpp
    Supporting material software to the article 
    STATISTICAL MECHANICS OF THE E COLI METABOLIC NETWORK IN STATIONARY GROWTH
    
    Copyright (C) 2016 D.De Martino, G.Tkacik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
                                                                           */


#include "iostream"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "fstream"
#include "vector"
#include "string"
using namespace std;

#define r 23                                // number of variables
#define N 86                               // number of constraints
#define SAMPLETIME 1e6 	                    // number of steps
#define samplerate 100                     // sampling rate
#define POLYHEDRON_FNAME "polcoli.dat"  // Name of the file containing the polytope                    
#define REACTIONS_BOUNDS_FNAME "boundscoli.dat"     // Name of the file containing the bounds of the reactions 
#define ELLIPSOID "ellis.dat"	// Name of the file containing the ellipsoid
#define INITPOINT "point0.dat"	     //  Name of the file containing the initial point
#define MAXFLUX 1e9			     // Initialization of variables tp, tm in extrema(). MAXFLUX should be larger than maximum flux
#define epsilon 1e-9
#define ID_RXN_TO_MAX 12                    //  flux maximization, reaction index, 12 is the biomass 
#define MAXIMUM_EXP 20
#define EPSILON	1e-10
#define BETA 0                           // flux maximization, rate of optimization, BETA=0: uniform sampling 
std::vector<int> matrice[N];                  // polytope matrix, pointers
std::vector<long double> matrice2[N];        // polytope matrix, stroichiometric coefficients
long double   boundmin[N];                // bounds
long double   boundmax[N]; 
long double   flux[r];                  // polytope variables
long double   versori[r][r];            // ellipsoid axis directions 
long double   cb[r];                 // optional: flux maximazation, objective function coefficients
long double   ran[2];               // random numbers
long double beta;	
long double previous[r];


void reading(){           // function that reads the input
    fstream file;
    string nomi;
    int rev;
    file.open(REACTIONS_BOUNDS_FNAME, ios::in);                                     // bounds reading
    for(int i=0;i<=N-1;i++){
                          file  >>  nomi  >> boundmin[i] >> boundmax[i];
                         
                           }   

    file.close();
    file.open(POLYHEDRON_FNAME, ios::in);                                     // polytope reading
    int  i1=0;
    int  j1=0;        
    long double pepe;
    while(file >> pepe){
	if(pepe!=0){
	    matrice[i1].push_back(j1);
	    matrice2[i1].push_back(pepe);
	}
	j1++;
	if(j1==r){
	    j1=0;
	    i1++;
	}
    }        
    file.close();  
   long double lenght;
   file.open(ELLIPSOID, ios::in);                       // ellipsoid reading
    for(int i=0;i<=r-1;i++){
            file >> lenght;
           
            for(int j=0;j<=r-1;j++) file >> versori[i][j]; 
            }
    file.close();
     file.open(INITPOINT, ios::in); 
    for(int i=0;i<=r-1;i++) file >> flux[i];
           
    file.close();
    for(int i=0;i<=r-1;i++){                  // obj function reading
                  cb[i]=0;
                  for(int j=0;j<=matrice[ID_RXN_TO_MAX].size()-1;j++)  cb[i] += versori[i][matrice[ID_RXN_TO_MAX][j]]*matrice2[ID_RXN_TO_MAX][j];
                  }
}


long double casual(){                  // a random number uniform in (0,1)
    long double   x = double  (random())/(RAND_MAX+1.0);
    return x;
}

void randomini(){           // random initialization
    int indip[r];
    int i1=0;
    for(int i=0;i<=N-1;i++){
	if(matrice[i].size()==1 && matrice2[i][0]==1 && i1<r){
	    indip[i1]=matrice[i][0];
	    i1++;
	}
    }
    for(int i=0;i<=r-1;i++) flux[i] = 0.01*boundmin[indip[i]]+0.01*casual()*(boundmax[indip[i]]-boundmin[indip[i]]);

}


void gaussMP(long double& a,long double& b){          // two gaussian random numbers
    long double   x,y,s;
    do{
	x=-1+2*casual();
	y=-1+2*casual();
	s = x*x+y*y;
    }while(s>=1);
    long double   c = sqrt(-2*log(s)/s);
    a = x*c;
    b = y*c;
}

void findextrema(long double& tp, long double& tm, int naxis){               // this function finds segment extrema along the direction of the axis "naxis"
    tp=MAXFLUX;
    tm=-MAXFLUX;
    for(int i=0;i<=N-1;i++){
       if(matrice[i].size()>0 && ( boundmin[i]!=-1000 || boundmax[i]!=1000)   ){
	long double x,y;
	x=y=0;
	for(int j=0;j<=matrice[i].size()-1;j++){
	    x += flux[matrice[i][j]]*matrice2[i][j];
	    y += versori[naxis][matrice[i][j]]*matrice2[i][j];
	}
       if(y!=0){
	long double   t1,t2;
	t1 = (boundmin[i]-x)/y;
	t2 = (boundmax[i]-x)/y; 
	if(t1<0 && t1>tm) tm=t1;
	if(t2<0 && t2>tm) tm=t2;
	if(t1>0 && t1<tp) tp=t1;
	if(t2>0 && t2<tp) tp=t2;

        }


      }
    }
}



void ellipsHR(){                    // hit and run with ellipsoid
    
    for(int i=0;i<=r-1;i++) previous[i]=flux[i];

    for(int yes=0;yes<=r-1;yes++){            // sweep over the axis directions
       
       long double   tp,tm;
       findextrema(tp,tm,yes);
       long double   t=0;
       int count=0;
       do{ 
         long double c = cb[yes];
         c*=BETA;
         long double ranvar = casual();
         if(fabs(c)<EPSILON) t = tm+(tp-tm)*ranvar;
         else{                                  //  flux maximization
	  

	   if(fabs(c*(tp-tm))<MAXIMUM_EXP) t = tm + log( 1 + ranvar*(exp(c*(tp-tm))-1))/c;
                                   
           else{
                    if(c<0) t = tm+log(1-ranvar)/c;
                    else    t = tp+log(ranvar)/c;
             }
        
           }
           count++;  
           }while((t-tm<EPSILON || tp-t<EPSILON) && count<10);

    if(count==10) t=0;
  
    for(int i=0;i<=r-1;i++) flux[i] += t*versori[yes][i];
   }
             int oko=1;
                       for(int i=0;i<=N-1;i++){
                                            double flux1=0;
                                            if(matrice[i].size()>0) for(int j=0;j<=matrice[i].size()-1;j++) flux1+=matrice2[i][j]*flux[matrice[i][j]];
                                           
                                           if(flux1<boundmin[i] || flux1>boundmax[i])oko=0;
                                            } 

                     if(oko==1)  for(int i=0;i<=r-1;i++) previous[i]=flux[i];
                     else for(int i=0;i<=r-1;i++) flux[i]=previous[i];

}


void minover(){   // optional: relaxation algorithm to find a point inside 
    int ok=0;
    long double alpha;
    long double xmin;

    do{

	int min;
	int sign;
	xmin=100000000000; 
	for(int i=0;i<=N-1;i++){
          if(matrice[i].size()>0){
	    long double x=0;
	    for(int j=0;j<=matrice[i].size()-1;j++) x += flux[matrice[i][j]]*matrice2[i][j];
	    long double x1 = x-boundmin[i];
	    long double x2 = boundmax[i]-x;
	    if(x1<xmin){
		xmin = x1;
		min =  i;
		sign = 1;
	    }
	    if(x2<xmin){
		xmin = x2;
		min =  i;
		sign = -1;
	    }
         }
	}
	if(xmin>0) ok=1;
	else{
	
           // cout << xmin << endl;
	    long double norm=0;
	    for(int j=0;j<=matrice[min].size()-1;j++) norm += matrice2[min][j]*matrice2[min][j];
	    alpha = -1.9*xmin/norm;

	    if(alpha<1e-16) alpha=1e-16;  
	    for(int j=0;j<=matrice[min].size()-1;j++) flux[matrice[min][j]] += sign*alpha*matrice2[min][j]; 
	}
    }while(ok==0);
   
  }







int main(){
    
    
    srand(time(0));         // seed for random numbers  
    reading();          // reading the input
    cout.precision(16);
    //randomini(); 
    minover();
   


    double previous[r];
    for(int i=0;i<=r-1;i++) previous[i]=flux[i];

      for(int i=0;i<=SAMPLETIME;i++){        // sampling   
	 ellipsHR();
         if(i%samplerate==0){
                                       for(int j=0;j<=N-1;j++){
                                              double value=0;  
                     
                                              for(int k=0;k<=matrice[j].size()-1;k++) value += matrice2[j][k]*flux[matrice[j][k]];
                                                         
                                                            
                                              cout << value << "  " ;
                                              }
                                       cout << endl;    
                    
                    
                      }
         }   
        
       
      
    return 0;
}
