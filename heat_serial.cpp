#include<iostream>
#include<stdlib.h>
#include<math.h> /* cos, pow */
#include<string.h> /* memcpy */
#include<fstream> /* ofstream */
using namespace std;
//#include "gnuplot-iostream.h"

int main(int argc, char *argv[]){

  const int size=atof(argv[1]);
  const double kappa = 1;
  //Array spacing is pi/length of array. deltay=deltax
  const double deltax=M_PI/size;
  //deltat should be smaller than deltax^2/(4*kappa)
  const double deltat=pow(deltax,2)/(10*kappa);
  const double time = 0.5*pow(M_PI,2)/kappa;
  //time step
  const int n = time/deltat;
  double temp_sum=0;
  double avg = 0;

  //Initializing square 2D array of pointers to array
  double** parray = new double*[size];
  for(int i = 0; i< size; i++)
    parray[i] = new double[size];

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
    parray[i][0] = pow(cos(i*deltax),2);
    parray[i][size-1] = pow(sin(i*deltax),2);
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
    for(int i = 1; i < size-1; i++){
      for(int j = 1; j < size-1; j++){
	parray_np[i][j] = parray[i][j]+deltat*kappa*(parray[i-1][j]+parray[i+1][j]+parray[i][j-1]+parray[i][j+1]-4*parray[i][j])/(pow(deltax,2));
      }
    }

    //Updating first and last row periodic boundary conditions
    for(int j=1; j < size-1; j++){
      parray_np[0][j]     =parray[0][j]+deltat*kappa*(parray[size-1][j]+parray[1][j]+parray[0][j-1]+parray[0][j+1]-4*parray[0][j])/(pow(deltax,2));
      parray_np[size-1][j]=parray[size-1][j]+deltat*kappa*(parray[size-2][j]+parray[0][j]+parray[size-1][j-1]+parray[size-1][j+1]-4*parray[size-1][j])/(pow(deltax,2));
    }

    //Moving array n+1 to n
    //memcpy(parray, parray_np, sizeof(parray));

    for(int i=0; i< size;i++){
      for(int j=0; j< size; j++){
	parray[i][j]=parray_np[i][j];
      }
    }
  }

  //Final volume averaged temperature
  for(int i = 0; i < size; i++){
      for(int j = 0; j < size; j++){
	temp_sum+=parray[i][j];
      }
  }

  avg=temp_sum/pow(size,2);

  cout<<"The volume averaged temperature is: " << avg << endl;

  //Writing final array to file
  ofstream dataFile;
  dataFile.open("temp_array.txt");
  dataFile << "Array Size: " << size <<"^2"<<endl;
  dataFile << "Volume averaged temperature: " << avg <<endl;
  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      dataFile << i <<" "<< j<<" "<< parray[i][j] <<endl;
    }
  }

  dataFile.close();
  
  cout<<"Success! Heat diffusion array written to file"<<endl;
  
  /* Gnuplot gp;
  gp << "set xrange [0:99]\n set yrange [0:99]\n";
  gp << "set pm3d\n";
  gp << "set hidden3d\n";
  gp << "set view map\n";
  gp << "splot '-' matrix" <<"\n";
  gp.send1d(parray);
  gp.flush();*/
}
