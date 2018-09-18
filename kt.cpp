#include <iostream>
#include <cmath>
#include <cstdlib>
#include<stdio.h>
#include<math.h>

#define X 1000
#define MR 0.125

#define L (3.14159265359*2)
#define T 1.2
#define dX ((double)L)/X

#define EX 1

#define RHO 1.0

using namespace std;

double f(double x);

double df(double u);

void KT(double u[], double dx, int nx);

void Ex1();

void Ex2();

void Ex3();

void printResults(double u[], char filename[]);

void plotGraph(double u[]);

int sign(double a) {
    if(a > 0)
        return 1;
    else if(a < 0)
        return -1;
    return 0;
}

double minVec(double vec[], int sizeVec) {
    double m = vec[0];
    for(int i=1; i<sizeVec; i++) {
        if(vec[i] < m)
            m = vec[i];
    }
    return m;
}

double maxVec(double vec[], int sizeVec) {
    double m = vec[0];
    for(int i=1; i<sizeVec; i++) {
        if(vec[i] > m)
            m = vec[i];
    }
    return m;
}

double minmod(double a, double b) {
    /*
    int i = 0;
    double vec[] = {a, b, c};
    while(i < 3) {
        if(vec[i] <= 0)
            break;
        i++;
    }
    if(i > 3)
        return minVec(vec, 3);
    i = 0;
    while(i < 3) {
        if(vec[i] >= 0)
            break;
        i++;
    }
    if(i > 3)
        return maxVec(vec, 3);
    return 0;
    */

    return 0.5*(sign(a) + sign(b)*min(abs(a), abs(b)));
}

void KT(double u[], double dx, int nx) {
    int J = nx-1;
    double ux[nx], uPlus[nx], uMinus[nx], a[nx], H1[nx], H2[nx], du[nx], solution[nx], dt;
    //double dx = (xJ-x0)/(nx-1);

    ///ux[0] = minmod(0, (u[1]-u[0])/(2.0*dx), theta*(u[1]-u[0])/(dx);
    /*ux[0] = 0;
    ux[J] = minmod(THETA * (u[J]-u[J-1])/dx,
                   (u[J+1]-u[J-1])/(2.0 * dx),
                   THETA * (u[J+1]-u[J])/dx);
    */

    ux[0] = minmod((u[0]-u[J])/dx, (u[1]-u[0])/dx);
    ux[J] = minmod((u[J]-u[J-1])/dx, (u[0]-u[J])/dx);

    for(int j = 1; j < nx; j++) {
        /*ux[j] = minmod(THETA*(u[j]-u[j-1])/dx,
                       (u[j+1]-u[j-1])/(2.0*dx),
                       THETA*(u[j+1]-u[j])/dx);
        */
        ux[j] = minmod((u[j]-u[j-1])/dx, (u[j+1]-u[j])/dx);
    }

    uMinus[0] = u[J] + double(dx/2.0)*ux[J];
    uPlus[0] = u[0] - double(dx/2.0)*ux[0];

    for(int j = 1; j < nx; j++) {
        uMinus[j] = u[j-1] + double(dx/2.0)*ux[j-1];
        uPlus[j] = u[j] - double(dx/2.0)*ux[j];
        a[j] = max(df(uMinus[j]), df(uPlus[j]));
    }

    for(int j = 0; j < J; j++) {
        H1[j] = double( (f(uPlus[j+1]) + f(uMinus[j+1]) - a[j+1]*(uPlus[j+1] - uMinus[j+1]))/2.0 );
    }
    H1[J] = double( (f(uPlus[0]) + f(uMinus[0]) - a[0]*(uPlus[0] - uMinus[0]))/2.0 );

    for(int j = 0; j < nx; j++) {
        H2[j] = double( (f(uPlus[j]) + f(uMinus[j]) - a[j]*(uPlus[j] - uMinus[j]))/2.0 );
    }

    for(int j = 0; j < nx; j++) {
        du[j] = (H2[j] - H1[j])/dx;
    }

    dt = MR*dx;
    for(int i = 0; i < nx; i++) {
        solution[i] = u[i] + dt*du[i];
    }

    plotGraph(solution);

}


void Ex1(){

  double u[X+3];
  int j;

  //initial conditions
  for(j = 1; j <= X+1; j++){
    u[j] = sin((j-1)*dX);
  }
  u[0] = u[1]; //ghost
  u[X+2] = u[X+1]; //ghost
  //-------------------

  KT(u, dX, X+3);
}

void Ex2(){

  double u[X+3];
  int j;

  //initial conditions
  for(j = 1; j <= X+1; j++){
    u[j] = sin((j-1)*dX) + 0.5;
  }
  u[0] = u[1]; //ghost
  u[X+2] = u[X+1]; //ghost
  //-------------------

  KT(u, dX, X+3);
}

void Ex3(){

  double u[X+3];
  int j;

  //initial conditions
  for(j = 1; j <= X+1; j++){
    u[j] = (j*dX < 1 ? 2 : -2);
  }
  u[0] = u[1]; //ghost
  u[X+2] = u[X+1]; //ghost
  //-------------------

  KT(u, dX, X+3);
}


int main()
{

    switch(EX){
    case 1: Ex1(); break;
    case 2: Ex2(); break;
    case 3: Ex3(); break;
    }

    return 0;
}

//########### f(u) ###########//
double f(double u){
  switch(EX){
    case 1: return u; break;
    case 2: return pow(u,2)/2; break;
    case 3: return (pow(u,2)-1)*(pow(u,2)-4)/4; break;
    default: return 0;
  }
}

double df(double u) {
    switch(EX){
    case 1: return 1; break;
    case 2: return u; break;
    case 3: return (pow(u,3) - double(2.5)*u); break;
    default: return 0;
  }
}

//###########################//
//           UTILS           //
//###########################//

void printResults(double u[], char filename[]){
  int j;
  FILE* f = fopen(filename, "w");
  for(j = 1; j <= X+1; j++){
    fprintf(f, "%lf\n", u[j]);
  }
  fclose(f);
}

void plotGraph(double u[]){
  FILE* temp = fopen("data.temp", "w");
  FILE* gnuplot = popen("gnuplot -persistent", "w");
  char* commands[] = {"set title \"TITLE\"", "plot 'data.temp' with lines"};
  int num_commands = 2;

  int j;
  for(j = 1; j <= X+1; j++){
    fprintf(temp, "%lf %lf \n", (j-1)*(dX), u[j]);
  }
  for(j = 0; j < num_commands; j++){
    fprintf(gnuplot, "%s \n", commands[j]);
  }
}

//###########################//
//           END ;)          //
//###########################//
