/*
    lovasz.cpp
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
#include "nrutil.h"
using namespace std;

#define r  23                                // number of variables
#define N  86                               // number of constraints
#define POLYHEDRON_FNAME "polcoli.dat"  // Name of the file containing the polyhedron                    
#define REACTIONS_BOUNDS_FNAME "boundscoli.dat"    // Name of the file containing the bounds of the reactions
#define INITPOINT "point0.dat"	     // 
#define MAXFLUX 10000.000000 			     // Initialization of variables tp, tm in extrema(). MAXFLUX should be larger than maximum flux
#define MINIMUM 1e-8

std::vector<int> matrice[N];
std::vector<long double> matrice2[N];
long double   boundmin[N];
long double   boundmax[N];
long double   flux[r]; 
long double   versori[r][r];
long double   lenght[r];

long double pythag(long double a,long double b){
    long double c =sqrt(a*a+b*b);
    return c;
}
void tred2(long double **a, int n, long double d[], long double e[])
    /*Householder reduction of a real, symmetric matrice a[1..n][1..n] . On output, a is replaced
      by the orthogonal matrice Q effecting the transformation. d[1..n] returns the diagonal ele-
      ments of the tridiagonal matrice, and e[1..n] the off-diagonal elements, with e[1]=0. Several
      statements, as noted in comments, can be omitted if only eigenvalues are to be found, in which
      case a contains no useful information on output. Otherwise they are to be included.*/
{
    int l,k,j,i;
    long double scale,hh,h,g,f;
    for (i=n;i>=2;i--) {
	l=i-1;
	h=scale=0.0;
	if (l > 1) {
	    for (k=1;k<=l;k++)
		scale += fabs(a[i][k]);
	    if (scale == 0.0){
		//Skip transformation.
		e[i]=a[i][l];
	    }
	    else {
		for (k=1;k<=l;k++) {
		    a[i][k] /= scale;
		    //Use scaled a’s for transformation.
		    h += a[i][k]*a[i][k];
		    //Form σ in h.
		}
		f=a[i][l];
		g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
		e[i]=scale*g;
		h -= f*g;
		//Now h is equation (11.2.4).
		a[i][l]=f-g;
		//Store u in the ith row of a.
		f=0.0;
		for (j=1;j<=l;j++) {
		    /* Next statement can be omitted if eigenvectors not wanted */
		    a[j][i]=a[i][j]/h;
		    //Store u/H in ith column of a.
		    g=0.0;
		    //Form an element of A · u in g.
		    for (k=1;k<=j;k++)
			g += a[j][k]*a[i][k];
		    for (k=j+1;k<=l;k++)
			g += a[k][j]*a[i][k];
		    e[j]=g/h;
		    //Form element of p in temporarily unused
		    //element of e.
		    f += e[j]*a[i][j];
		}
		hh=f/(h+h);
		//Form K, equation (11.2.11).
		for (j=1;j<=l;j++) {
		    //Form q and store in e overwriting p.
		    f=a[i][j];
		    e[j]=g=e[j]-hh*f;
		    for (k=1;k<=j;k++)
			//Reduce a, equation (11.2.13).
			a[j][k] -= (f*e[k]+g*a[i][k]);
		}
	    }
	} else
	    e[i]=a[i][l];
	d[i]=h;
    }
    /* Next statement can be omitted if eigenvectors not wanted */
    d[1]=0.0;
    e[1]=0.0;
    /* Contents of this loop can be omitted if eigenvectors not
       wanted except for statement d[i]=a[i][i]; */
    for (i=1;i<=n;i++) {
	//Begin accumulation of transformation matrices
	l=i-1;
	if (d[i]) {
	    //This block skipped when i=1.
	    for (j=1;j<=l;j++) {
		g=0.0;
		for (k=1;k<=l;k++)
		    //Use u and u/H stored in a to form P·Q.
		    g += a[i][k]*a[k][j];
		for (k=1;k<=l;k++)
		    a[k][j] -= g*a[k][i];
	    }
	}
	d[i]=a[i][i];
	//This statement remains.
	a[i][i]=1.0;
	//Reset row and column of a to identity
	for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	//matrice for next iteration.
    }
}

void tqli(long double d[], long double e[], int n, long double **z)
    /*QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, sym-
      metric, tridiagonal matrice, or of a real, symmetric matrice previously reduced by tred2 §11.2. On
      input, d[1..n] contains the diagonal elements of the tridiagonal matrice. On output, it returns
      the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrice,
      with e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues, several lines
      may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrice are de-
      sired, the matrice z[1..n][1..n] is input as the identity matrice. If the eigenvectors of a matrice
      that has been reduced by tred2 are required, then z is input as the matrice output by tred2.
      In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].*/
{
    long double pythag(long double a, long double b);
    int m,l,iter,i,k;
    long double s,r1,p,g,f,dd,c,b;
    for (i=2;i<=n;i++) e[i-1]=e[i];
    //Convenient to renumber the el-
    e[n]=0.0;
    //ements of e.
    for (l=1;l<=n;l++) {
	iter=0;
	do {
	    for (m=l;m<=n-1;m++) {
		//Look for a single small subdi-
		dd=fabs(d[m])+fabs(d[m+1]);
		//agonal element to split
		if ((long double)(fabs(e[m])+dd) == dd) break;
		//the matrice.
	    }
	    if (m != l) {
		//if (iter++ == 30) nrerror("Too many iterations in tqli");
		g=(d[l+1]-d[l])/(2.0*e[l]);
		//Form shift.
		r1=pythag(g,1.0);
		g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
		//This is dm − ks .
		s=c=1.0;
		p=0.0;
		for (i=m-1;i>=l;i--) {
		    //A plane rotation as in the origi-
		    f=s*e[i];
		    //nal QL, followed by Givens
		    b=c*e[i];
		    //rotations to restore tridiag-
		    e[i+1]=(r1=pythag(f,g));
		    //onal form.
		    if (r1 == 0.0) {
			//Recover from underflow.
			d[i+1] -= p;
			e[m]=0.0;
			break;
		    }
		    s=f/r1;
		    c=g/r1;
		    g=d[i+1]-p;
		    r1=(d[i]-g)*s+2.0*c*b;
		    d[i+1]=g+(p=s*r1);
		    g=c*r1-b;
		    /* Next loop can be omitted if eigenvectors not wanted*/
		    for (k=1;k<=n;k++) {
			//Form eigenvectors.
			f=z[k][i+1];
			z[k][i+1]=s*z[k][i]+c*f;
			z[k][i]=c*z[k][i]-s*f;
		    }
		}
		if (r1 == 0.0 && i >= l) continue;
		d[l] -= p;
		e[l]=g;
		e[m]=0.0;
	    }
	} while (m != l);
    }
}

void reading(){
    fstream file;
    string nomi;
    file.open(REACTIONS_BOUNDS_FNAME, ios::in);
    int rev ;                                    // bounds reading
    for(int i=0;i<=N-1;i++){
                       file >> nomi >> boundmin[i] >> boundmax[i];
                       
                       }     
                           
    file.close();
    file.open(POLYHEDRON_FNAME, ios::in);                                     // polyhedron reading
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

    // file.open(INITPOINT, ios::in); 
       for(int i=0;i<=r-1;i++) flux[i]=0;
           
       
     //  for(int i=0;i<=r-1;i++) flux[i]=0;

      
  
}






void inizialization(){

       
    //  for(int i=0;i<=r-1;i++) flux[i] = 0;
     


     /* minover();
      cout << "found point inside" << endl;*/
      for(int i=0;i<=r-1;i++) for(int j=0;j<=r-1;j++){
                                                             if(i==j) versori[i][j]=1;
                                                             else     versori[i][j]=0;
                                                      }
      long double maxi = 0;
      for(int i=0;i<=N-1;i++) if(boundmax[i]-boundmin[i]>maxi) maxi = boundmax[i]-boundmin[i];
      long double Radius =sqrt(r)*maxi;
     // Radius=10000;
      for(int i=0;i<=r-1;i++) lenght[i] = Radius;



}

void calcvolume(){
        long double V=0;
        for(int i=0;i<=r-1;i++) V+=log(lenght[i]);
      //  cout << " logvolume  " << V << endl;
        }
void basicellipsoid(int value){
    
                    
   
   
    
    long double** Matr = new long double*[r+1];
    for(int i = 0; i < r+1; ++i) Matr[i] = new long double[r+1];
    
    for(int i=0;i<=r;i++) for(int j=0;j<=r;j++) Matr[i][j]=0;
 
    for(int k=0;k<=r-1;k++) for(int i=0;i<=r-1;i++) for(int j=0;j<=r-1;j++) Matr[i+1][j+1] += versori[k][i]*versori[k][j]*lenght[k]*lenght[k]; 


    long double fluxstep[r];
    long double norm=0;
    for(int i=0;i<=r-1;i++) fluxstep[i]=0;
    int sign=1;
    if(value>=N) sign=-1;    

    if(sign==1){  
                  long double newvec[r];
                  for(int i=0;i<=r-1;i++) newvec[i]=0;
                  for(int i=0;i<=matrice[value].size()-1;i++) newvec[matrice[value][i]]=-matrice2[value][i];
                  for(int i=0;i<=r-1;i++){
                                  long double x=0;
                                  for(int j=0;j<=r-1;j++) x += Matr[i+1][j+1]*newvec[j];
                                  fluxstep[i]=x;
                                  norm += x*newvec[i];
                                  }         
                }
    else{ 
                  long double newvec[r];
                  for(int i=0;i<=r-1;i++) newvec[i]=0;
                  for(int i=0;i<=matrice[value-N].size()-1;i++) newvec[matrice[value-N][i]]=matrice2[value-N][i];
                  for(int i=0;i<=r-1;i++){
                                  long double x=0;
                                  for(int j=0;j<=r-1;j++) x += Matr[i+1][j+1]*newvec[j];
                                  fluxstep[i]=x;
                                  norm += x*newvec[i];
                                  }    
                  



                 }
   
    norm=sqrt(norm);
    for(int i=0;i<=r-1;i++) fluxstep[i]/=norm;
    long double rr=r;
    long double a1,a2,a3;
    a1 = 1./(rr+1.);
    
    for(int i=0;i<=r-1;i++) flux[i]   -=   a1*fluxstep[i];
                           
                        
    a2 = rr*rr/(rr*rr-1.);
    a3 = 2.;
       
   



    
    



    for(int i=1;i<=r;i++) for(int j=1;j<=r;j++) Matr[i][j] =    a2*(Matr[i][j]-a3*a1*fluxstep[i-1]*fluxstep[j-1]);   // assign matrix shallow cut

                                       
                                       
                                       
                                        
    long double D[r+1];                                                                             
    long double E[r+1];
    tred2(Matr,r,D,E); 
    tqli(D,E,r,Matr);
    for(int i=0;i<=r-1;i++){    
                               // cout << lenght[i] <<"  ";
                             //   lenght[i] = sqrt(D[i+1]);
                                                // attenzio, piccola correzio!
                                  if(D[i+1]<0) D[i+1]=-D[i+1];
                                  lenght[i]=sqrt(D[i+1]);
                               // cout << lenght[i] << endl;
                                
                          }
    for(int i=1;i<=r;i++) for(int j=1;j<=r;j++) versori[i-1][j-1]=Matr[j][i];
    for(int i=0;i<=r;i++) free(Matr[i]);
    free(Matr);
}

void newellipsoid(int value){
    
                    
   long double beta;
   int sign=1;
   if(value>=N) sign=-1;
   if(sign==1){
           beta=-boundmin[value];
           for(int i=0;i<=matrice[value].size()-1;i++) beta += matrice2[value][i]*flux[matrice[value][i]];
           }
   else{
        beta=boundmax[value-N];
        for(int i=0;i<=matrice[value-N].size()-1;i++) beta -= matrice2[value-N][i]*flux[matrice[value-N][i]];
        }
   
    
    long double** Matr = new long double*[r+1];
    for(int i = 0; i < r+1; ++i) Matr[i] = new long double[r+1];
    
    for(int i=0;i<=r;i++) for(int j=0;j<=r;j++) Matr[i][j]=0;
 
    for(int k=0;k<=r-1;k++) for(int i=0;i<=r-1;i++) for(int j=0;j<=r-1;j++) Matr[i+1][j+1] += versori[k][i]*versori[k][j]*lenght[k]*lenght[k]; 


    long double fluxstep[r];
    long double norm=0;
    for(int i=0;i<=r-1;i++) fluxstep[i]=0;
    
    if(sign==1){  
                  long double newvec[r];
                  for(int i=0;i<=r-1;i++) newvec[i]=0;
                  for(int i=0;i<=matrice[value].size()-1;i++) newvec[matrice[value][i]]=-matrice2[value][i];
                  for(int i=0;i<=r-1;i++){
                                  long double x=0;
                                  for(int j=0;j<=r-1;j++) x += Matr[i+1][j+1]*newvec[j];
                                  fluxstep[i]=x;
                                  norm += x*newvec[i];
                                  }         
                }
    else{ 
                  long double newvec[r];
                  for(int i=0;i<=r-1;i++) newvec[i]=0;
                  for(int i=0;i<=matrice[value-N].size()-1;i++) newvec[matrice[value-N][i]]=matrice2[value-N][i];
                  for(int i=0;i<=r-1;i++){
                                  long double x=0;
                                  for(int j=0;j<=r-1;j++) x += Matr[i+1][j+1]*newvec[j];
                                  fluxstep[i]=x;
                                  norm += x*newvec[i];
                                  }    
                  



                 }
   
    norm=sqrt(norm);
    beta/=norm;
    long double rr=r;
    long double a1,a2,a3;
    a1 = (1.-rr*beta)/(rr+1.);
   
    for(int i=0;i<=r-1;i++){
                     
                     fluxstep[i]/=norm;
                     flux[i]   -=   a1*fluxstep[i];
                     
                       }
    a2 = rr*rr*(1.-beta*beta)/(rr*rr-1.);
    a3 = 2./(1.-beta);
       
  //  cout << " beta " << beta << "  " << a1 << "   " << norm << endl;



    
    



    for(int i=1;i<=r;i++) for(int j=1;j<=r;j++){

                                       // cout << i << "  " << j << "  " << Matr[i][j]<< "  ";
                                        Matr[i][j] =    a2*(Matr[i][j]-a3*a1*fluxstep[i-1]*fluxstep[j-1]);   // assign matrix shallow cut
                                        //cout << Matr[i][j] << endl;
                                        }
    long double D[r+1];                                                                             
    long double E[r+1];
    tred2(Matr,r,D,E); 
    tqli(D,E,r,Matr);
    for(int i=0;i<=r-1;i++){    
                               // cout << lenght[i] <<"  ";
                                                // attenzio, piccola correzio!
                                  if(D[i+1]<0) D[i+1]=-D[i+1];
                                lenght[i] = sqrt(D[i+1]);

                                 // lenght[i]=D[i+1];
                               // cout << lenght[i] << endl;
                                
                          }
    for(int i=1;i<=r;i++) for(int j=1;j<=r;j++) versori[i-1][j-1]=Matr[j][i];
    for(int i=0;i<=r;i++) free(Matr[i]);
    free(Matr);
}

int checkinside(){
                   int check=-1;
                   int r1=0;
                   do{         
                       
                      int min,sign;
                      long double xmin=10000;
                      int i=0;  
                      do{
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
                                       i++;
	                               }while(i<N && xmin>=0);
                       if(xmin<0){
                                   check =  min;
                                   if(sign<0) check+=N;   
                                  // cout << xmin << endl;    
                                 }
                       else r1++;
                       
                       }while(r1<r && check<0);        
                     

                       return check;

}



int checkaxis(){
                   int check=-1;
                   int r1=0;   
                   long double div = r+1;
                   do{         
                      long double punto1[r];
                      long double punto2[r];
                      for(int i=0;i<=r-1;i++){
                                              punto1[i]=punto2[i]=flux[i];
                                              long double step = lenght[r1]*versori[r1][i]/div;
                                              punto1[i]+=step;
                                              punto2[i]-=step;
                                             }   
                      int min,sign;
                      long double xmin=10000;  
                      int i=0;
                      do{
                                      if(matrice[i].size()>0){
	                                               long double x=0;
	                                               for(int j=0;j<=matrice[i].size()-1;j++) x += punto1[matrice[i][j]]*matrice2[i][j];
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
                                                       x=0;
                                                       for(int j=0;j<=matrice[i].size()-1;j++) x += punto2[matrice[i][j]]*matrice2[i][j];
                                                       x1 = x-boundmin[i];
	                                               x2 = boundmax[i]-x;
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
                                       i++;
	                               }while(i<N && xmin>=0);
                       if(xmin<0){
                                   check =  min;
                                   if(sign<0) check += N;   
                                  // cout << xmin << endl;    
                                 }
                       else r1++;
                       
                       }while(r1<r && check<0);        
                       

                       return check;

}




int main(){


            reading();
            inizialization();
            int check; 
            long double tstart, tstop, ttime;
            tstart = (long double)clock()/CLOCKS_PER_SEC;
    
         
   
            do{
                   check=1000;
                   int check0 = checkinside();
                   if(check0<0){
                          //       cout << "inside" << endl;
                                 check =  checkaxis();
                                 if(check>=0) newellipsoid(check);
                                  calcvolume();
                                 }
                   else{
                              basicellipsoid(check0);
                            //  cout << "outside" << endl;
                             
       
                      }
              }while(check>=0);
              tstop = (long double)clock()/CLOCKS_PER_SEC;
             ttime= tstop-tstart; 
            // cout << "result" << endl;
             for(int i=0;i<=r-1;i++) cout << flux[i] << "  ";
             cout << endl;
             for(int i=0;i<=r-1;i++){
                        cout << lenght[i] << endl;
                          for(int j=0;j<=r-1;j++) cout << versori[i][j] << "  ";
                         cout << endl;
                         }           
              //  cout << ttime << endl;
    return 0;
}
