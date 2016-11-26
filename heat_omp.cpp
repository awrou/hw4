#include<iostream>
#include<stdlib.h>
#include<math.h> /* cos, pow */
#include<string.h> /* memcpy */
#include<fstream> /* ofstream */
#include<omp.h>
using namespace std;
//#include "gnuplot-iostream.h"

double sqr(double x){
  return x*x;
}

int main(int argc, char *argv[]){
  
  double t1 = omp_get_wtime();
  
  const int size = atof(argv[1]);
  const int nthreads = atof(argv[2]);
  
  const double kappa = 1;
  //Array spacing is pi/length of array. deltay=deltax
  const double deltax=M_PI/size;
  //deltat should be smaller than deltax^2/(4*kappa)
  const double deltat = sqr(deltax)/(6*kappa);
  const double time = 0.5*sqr(M_PI)/kappa;
  //time step
  const int n = time/deltat;
  double temp_sum=0;
  double avg = 0;

  //Initializing square 2D array of pointers to array
  double** parray = new double*[size];
  for(int i = 0; i< size; i++){
    parray[i] = new double[size];
  }
  //Initializing n+1 square 2D array of pointers to array
  double** parray_np = new double*[size];
  for(int i = 0; i< size; i++){
    parray_np[i] = new double[size];
  }
  //Zeroing array
  for(int i = 0; i< size; i++){
    for(int j = 0; j< size; j++){
      parray[i][j]=0;
    }
  }
  //Initializing array boundary conditions
  for(int i = 0; i< size; i++){
    //First and last column initialized
    parray[i][0] = sqr(cos(i*deltax));
    parray[i][size-1] = sqr(sin(i*deltax));
    //parray[0][i] = parray[size-1][i];
  }

  //Equating arrays at t0
  // memcpy(parray_np, parray, sizeof(parray));
  for(int i=0; i< size;i++){
    for(int j=0; j< size; j++){
      parray_np[i][j]=parray[i][j];
    }
  }

  //Doing n time steps of updates
  for(int k = 0; k < n; k++){
    //Updating center blocks
    #pragma omp parallel for num_threads(nthreads)
    for(int i = 1; i < size-1; i++){
      for(int j = 1; j < size-1; j++){
	parray_np[i][j] = parray[i][j]+deltat*kappa*(parray[i-1][j]+parray[i+1][j]+parray[i][j-1]+parray[i][j+1]-4*parray[i][j])/(sqr(deltax));
      }
    }
    
    //Updating first and last row periodic boundary conditions
    #pragma omp parallel for num_threads(nthreads)
    for(int j=1; j < size-1; j++){
      parray_np[0][j]     =parray[0][j]+deltat*kappa*(parray[size-1][j]+parray[1][j]+parray[0][j-1]+parray[0][j+1]-4*parray[0][j])/(sqr(deltax));
      parray_np[size-1][j]=parray[0][j];
    }

    //Moving array n+1 to n
    //memcpy(parray, parray_np, sizeof(parray));
    #pragma omp parallel for num_threads(nthreads)
    for(int i=0; i< size;i++){
      for(int j=0; j< size; j++){
	parray[i][j]=parray_np[i][j];
      }
    }
  }

  //Final volume averaged temperature
  //#pragma omp parallel for num_threads(nthreads)
  for(int i = 0; i < size; i++){
      for(int j = 0; j < size; j++){
	temp_sum+=parray[i][j];
      }
  }

  avg=temp_sum/sqr(size);

  double t2 = omp_get_wtime() -t1;
  
  cout<<"The volume averaged temperature is: " << avg << endl;
  cout<<"Total Runtime: "<< t2 <<" s" <<endl;

  //Writing final array to file
  ofstream dataFile;
  char buffer[33];
  sprintf(buffer,"temp_array_omp_%d.txt",size);
  dataFile.open(buffer);
  dataFile << "Array Size: " << size <<"^2"<<endl;
  dataFile << "Volume averaged temperature: " << avg <<endl;
  dataFile << "Total Runtime: " << t2 <<" s" <<endl;
  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      dataFile << i <<" "<< j<<" "<< parray[i][j] <<endl;
    }
  }

  dataFile.close();
  
  cout<<"Success! Heat diffusion data written to file"<<endl;
  
  /* Gnuplot gp;
  gp << "set xrange [0:99]\n set yrange [0:99]\n";
  gp << "set pm3d\n";
  gp << "set hidden3d\n";
  gp << "set view map\n";
  gp << "splot '-' matrix" <<"\n";
  gp.send1d(parray);
  gp.flush();*/
}
