#include<stdio.h>
#include<math.h>

#define N 400
#define EX 1

void firstorder();

void thirdorder();

double calcUf(double U[N], int i);

void init(double u[N], double* tf, double dx);

void printResults(double u[N], char filename[], int t);

int main(){
  firstorder();
  thirdorder();
  return 0;
}

void firstorder(){
  double u[N], up[N];
  double a = 1;
  double dt = 0.00025; //0.00025
  double dx = 0.005;
  double x, t, tf;
  int i, j;

  init(u, &tf, dx);

  //main loop
  for(i = 1, t = dt; t < 1; i++){
    for(j = 0; j < N; j++){
      up[j] = u[j];
    }
    t = i*dt;
    for(j = 1; j < N-1; j++){
      x = j*dx;
      u[j] = up[j]-(dt*a/dx)*(up[j]-up[j-1]);
    }
    printResults(u, "resultsFOU.txt", i);
  }
  printf("Finish at t = %.2f\n", t-dt);
}

void thirdorder(){
  double u[N], up[N];
  double a = 1;
  double dt = 0.00025; //0.00025
  double dx = 0.005;
  double x, t, tf;
  double uf, ug;
  double aplus, aminus;
  int i, j;

  init(u, &tf, dx);

  //main loop
  for(i = 1, t = dt; t < 1; i++){
    for(j = 0; j < N; j++){
      up[j] = u[j];
    }
    t = i*dt;
    //using FOU to get u[1]
    u[1] = up[1]-(dt*a/dx)*(up[1]-up[0]);
    uf = calcUf(up, 1);
    for(j = 2; j < N-1; j++){
      //for each x
      x = j*dx;
      ug = uf;
      uf = calcUf(up, j);
      u[j] = up[j]-(dt*a/dx)*(uf-ug);
    }
    printResults(u, "resultsTOPUS.txt", i);
  }
  printf("Finish at t = %.2f\n", t-dt);
}

double calcUf(double U[N], int i){
  double ufhat, uuhat;
  if(U[i+1]-U[i-1] != 0){
    uuhat = (U[i]-U[i-1])/(U[i+1]-U[i-1]);
    if(0 <= uuhat && uuhat <= 1){
      ufhat = 2*pow(uuhat, 4)-3*pow(uuhat, 3)+2*uuhat;
    } else {
      ufhat = uuhat;
    }
    return (ufhat*(U[i+1]-U[i-1]) + U[i-1]);
  } else {
    return U[i];
  }
}

void printResults(double u[N], char filename[], int t){
  int k;
  FILE* f = fopen(filename, "a");
  fprintf(f, "%d ", t);
  for(k = 0; k < N; k++){
    fprintf(f, "%lf ", u[k]);
  }
  fprintf(f, "\n");
  fclose(f);
}

void init(double u[N], double* tf, double dx){
  int i;
  double x;

  //initial conditions
  if(EX == 1){
    for(i = 0; i < N; i++){
      x = i*dx;
      if(x < 0.2){
        u[i] = exp(-log(50)*pow((x-0.15)/0.05, 2));
      } else if(0.3 < x && x < 0.4){
        u[i] = 1;
      } else if(0.5 < x && x < 0.55){
        u[i] = 20*x-10;
      } else if(0.55 <= x && x < 0.66){
        u[i] = -20*x+12;
      } else if(0.7 < x && x < 0.8){
        u[i] = sqrt(1-pow((x-0.75)/0.05, 2));
      } else {
        u[i] = 0;
      }
    }
    *tf = 1;
  } else {
    for(i = 0; i < N; i++){
      x = i*dx;
      if(x < 1){
        u[i] = 0;
      } else if(1 <= x && x <= 1.2){
        u[i] = 1;
      } else if(1.2 < x && x <= 1.4){
        u[i] = 4*x-0.6;
      } else if(1.4 < x && x <= 1.6){
        u[i] = (-4)*x+2.6;
      } else if(1.6 < x && x <= 1.8){
        u[i] = 1;
      } else {
        u[i] = 0;
      }
    }
    *tf = 0.25;
  }

  //boundary conditions
  u[0] = u[N-1] = 0;
}

void printVTK(double u[N], int t){
  /*int j, k;
  char filename[30];
  sprintf(filename, "tempFinal/e%d.vtk", t);
  FILE* f = fopen(filename, "w");
  fprintf(f, "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS 64 64 1\nX_COORDINATES 64 double\n0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 \nY_COORDINATES 64 double\n0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 \nZ_COORDINATES 1 double\n0\nPOINT_DATA 4096\nFIELD FieldData 1\nsolucao 1 4096 double\n");
  for(k = 1; k <= 64; k++){
    for(j = 1; j <= 64; j++){
      fprintf(f, "%lf ", u[j][k]);
    }
  }
  fprintf(f, "\n");
  fclose(f);*/
}
