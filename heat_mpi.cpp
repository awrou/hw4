#include<iostream>
#include<stdlib.h>
#include<math.h> /* cos, pow */
#include<string.h> /* memcpy */
#include<fstream> /* ofstream */
#include<time.h> /* clock_t, CLOCKS_PER_SEC */
#include<mpi.h>
using namespace std;
//#include "gnuplot-iostream.h"

double sqr(double x){
  return x*x;
}

int main(int argc, char *argv[]){
  
  clock_t t;
  t = clock();
  const int size=atof(argv[1]);
  const double kappa = 1;
  //Array spacing is pi/length of array. deltay=deltax
  const double deltax=M_PI/size;
  //deltat should be smaller than deltax^2/(4*kappa)
  const double deltat=sqr(deltax)/(6*kappa);
  const double time = 0.5*sqr(M_PI)/kappa;
  //time step
  const int n = time/deltat;
  double temp_sum=0;
  double avg = 0;

  double* halofs=new double[size];
  double* halobs=new double[size];
  double* halofr=new double[size];
  double* halobr=new double[size];

  //MPI Init
  int nprocess;
  int rank;
  MPI_Status status[4];
  MPI_Request request[4];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nrows = size/nprocess;

  //Initializing square 2D array of pointers to array
  double** parray = new double*[nrows + 2];
  for(int i = 0; i< size; i++){
    parray[i] = new double[size];
  }
  
  //Initializing n+1 square 2D array of pointers to array
  double** parray_np = new double*[nrows + 2];
  for(int i = 0; i< nrows+2; i++){
    parray_np[i] = new double[size];
  }
  //Zeroing array
  for(int i = 0; i< nrows+2; i++){
    for(int j = 0; j< size; j++){
      parray[i][j]=0;
    }
  }
  //Initializing array boundary conditions
  for(int i = 0; i< nrows; i++){
    //First and last column initialized
    parray[i+1][0] = sqr(cos((i+rank*nrows)*deltax));
    parray[i+1][size-1] = sqr(sin((i+rank*nrows)*deltax));
    //parray[0][i] = parray[size-1][i];
  }

  //Equating arrays at t0
  // memcpy(parray_np, parray, sizeof(parray));
  for(int i=0; i< nrows+2;i++){
    for(int j=0; j< size; j++){
      parray_np[i][j]=parray[i][j];
    }
  }

  //Doing n time steps of updates
  for(int k = 0; k < n; k++){
    //Updating center blocks
    for(int i = 1; i < nrows+1; i++){
      for(int j = 1; j < size-1; j++){
	parray_np[i][j] = parray[i][j]+deltat*kappa*(parray[i-1][j]+parray[i+1][j]+parray[i][j-1]+parray[i][j+1]-4*parray[i][j])/(sqr(deltax));
      }
    }

    //Moving array n+1 to n
    //memcpy(parray, parray_np, sizeof(parray));
    for(int i=0; i< nrows+2;i++){
      for(int j=0; j< size; j++){
	parray[i][j]=parray_np[i][j];
      }
    }
    for(int i=0; i< size; i++){
      halobs[i] = parray[1][i];
      halofs[i] = parray[nrows][i];
    }

    int rf=(rank+1)%nprocess;
    int rb=(rank-1)%nprocess;
    if(rb<0){
      rb=nprocess-1;
    }

    MPI_Isend(&halofs, size, MPI_DOUBLE, rf, 1, MPI_COMM_WORLD,&request[0]);
    MPI_Isend(&halobs, size, MPI_DOUBLE, rb, 2, MPI_COMM_WORLD,&request[1]);
    MPI_Irecv(&halofr, size, MPI_DOUBLE, rf, 2, MPI_COMM_WORLD,&request[2]);
    MPI_Irecv(&halobr, size, MPI_DOUBLE, rb, 1, MPI_COMM_WORLD,&request[3]);

    MPI_Waitall(4, request, status);

    for(int i=0;i< size;i++){
      parray[0][i]=halofr[i];
      parray[nrows+1][i]=halobr[i];
    }

    //Final volume averaged temperature
    for(int i = 0; i < size; i++){
      for(int j = 0; j < size; j++){
	temp_sum+=parray[i][j];
      }
    }
  }

  avg=temp_sum/sqr(size);

  t =(clock() - t);
  
  cout<<"The volume averaged temperature is: " << avg << endl;
  cout<<"Total Runtime: "<<double(t)/CLOCKS_PER_SEC <<" s" <<endl;

  //Writing final array to file
  ofstream dataFile;
  char buffer[33];
  sprintf(buffer,"temp_array_mpi_%d.txt",size);
  dataFile.open(buffer);
  dataFile << "Array Size: " << size <<"^2"<<endl;
  dataFile << "Volume averaged temperature: " << avg <<endl;
  dataFile << "Total Runtime: " << double(t)/CLOCKS_PER_SEC <<" s" <<endl;
  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      dataFile <<i<<" "<<j<<" "<< parray[i][j] <<endl;
    }
  }

  dataFile.close();
  
  cout<<"Success! Heat diffusion data written to file"<<endl;

  MPI_Finalize();

  /* Gnuplot gp;
  gp << "set xrange [0:99]\n set yrange [0:99]\n";
  gp << "set pm3d\n";
  gp << "set hidden3d\n";
  gp << "set view map\n";
  gp << "splot '-' matrix" <<"\n";
  gp.send1d(parray);
  gp.flush();*/
  
}
