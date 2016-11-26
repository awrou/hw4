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

  //MPI Init
  int nprocess;
  int rank;
  MPI_Status status[4];
  MPI_Request request[4];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  double t1, t2;
  t1 = MPI_Wtime();
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

  int nrows = size/nprocess;

  //Initializing square 2D array of pointers to array
  double parray[nrows+2][size];
    for(int i=0;i<nrows+2;i++){
      for(int j=0; j<size;j++){
	parray[i][j]=0;
      }
    }

  //Initializing array boundary conditions
    for(int i = 0; i< size; i++){
	parray[1][i] = sqr(cos((i+rank*nrows)*deltax));
	parray[size/nprocess][i] = sqr(sin((i+rank*nrows)*deltax));
      }

  double parray_np[nrows+2][size];
    for(int i=0;i<nrows+2;i++){
      for(int j=0; j<size;j++){
	parray_np[i][j]=parray[i][j];
      }
    }

  //Doing n time steps of updates
  for(int k = 0; k < n; k++){
    int rowf =1;
    if(rank==0){ rowf=2; } 
    int rowe = nrows+1;
    if(rank==nprocess-1){ rowe = nrows; }

    //Updating center blocks
    for(int i = rowf; i < rowe; i++){
      for(int j = 1; j < size-1; j++){
	parray_np[i][j] = parray[i][j]+deltat*kappa*(parray[i-1][j]+parray[i+1][j]+parray[i][j-1]+parray[i][j+1]-4*parray[i][j])/(sqr(deltax));
      }
    }

    for(int i=rowf; i < rowe; i++){
      parray_np[i][size-1]     =parray[i][size-1]+deltat*kappa*(parray[i][size-2]+parray[i][0]+parray[i-1][size-1]+parray[i+1][size-1]-4*parray[i][size-1])/(sqr(deltax));
      parray_np[i][0]=parray_np[i][size-1];
    }

    //Moving array n+1 to n
    for(int i=0; i< nrows+2;i++){
      for(int j=0; j< size; j++){
	parray[i][j]=parray_np[i][j];
      }
    }

    int rf=(rank+1) % nprocess;
    int rb=(rank-1) % nprocess;
    if(rb<0){
      rb=nprocess +rb;
    }

    MPI_Isend(&parray[size/nprocess][0], size, MPI_DOUBLE, rf, 1, MPI_COMM_WORLD,&request[0]);
    MPI_Isend(&parray[1][0], size, MPI_DOUBLE, rb, 2, MPI_COMM_WORLD,&request[1]);
    MPI_Irecv(&parray[size/nprocess+1][0], size, MPI_DOUBLE, rf, 2, MPI_COMM_WORLD,&request[2]);
    MPI_Irecv(&parray[0][0], size, MPI_DOUBLE, rb, 1, MPI_COMM_WORLD,&request[3]);

    MPI_Waitall(4, request, status);

    //Final volume averaged temperature
    for(int i = 1; i < nrows+1; i++){
      for(int j = 0; j < size; j++){
   	temp_sum+=parray[i][j];
      }
    }
    
  }

  avg=temp_sum/sqr(size);

  t2 = MPI_Wtime() - t1;

  MPI_Finalize();
  
  cout<<"The volume averaged temperature is: " << avg << endl;
  cout<<"Total Runtime: "<< t2 <<" s" <<endl;

  //Writing final array to file
  ofstream dataFile;
  char buffer[33];
  sprintf(buffer,"temp_array_mpi_%d.txt",size);
  dataFile.open(buffer);
  dataFile << "Array Size: " << size <<"^2"<<endl;
  dataFile << "Volume averaged temperature: " << avg <<endl;
  dataFile << "Total Runtime: " << t2 <<" s" <<endl;
  for(int i=1; i<nrows+1; i++){
    for(int j=0; j<size; j++){
      dataFile <<i<<" "<<j<<" "<< parray[i][j] <<endl;
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
