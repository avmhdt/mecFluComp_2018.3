#include<stdio.h>
#include<math.h>

#define X 1000
#define MR 0.125

#define L 2//(3.14159265359*2)
#define T 1.2
#define dX ((double)L)/X

#define EX 3

double f(double x);

void LxF(double u[]);

void Ex1();

void Ex2();

void Ex3();

void printResults(double u[], char filename[]);

void plotGraph(double u[]);

int main(){
  switch(EX){
    case 1: Ex1(); break;
    case 2: Ex2(); break;
    case 3: Ex3(); break;
  }
  return 0;
}

void LxF(double u[]){
  double up[X+3];
  double dT = MR*(dX);
  double lbd = MR;
  int t, j;

  for(t = 1; t*dT <= T; t++){
    for(j = 1; j <= X+1; j++){
      up[j] = u[j];
    }
    up[0] = u[1];
    up[X+2] = u[X+1];

    for(j = 1; j <= X+1; j++){
      u[j] = (up[j+1]+up[j-1])/2.0 - (lbd/2.0)*(f(up[j+1])-f(up[j-1]));
    }
  }
  //printResults(u, "Ex1LxFhmm.txt");
  plotGraph(u);
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

  LxF(u);
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

  LxF(u);
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

  LxF(u);
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
//###########################//


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
