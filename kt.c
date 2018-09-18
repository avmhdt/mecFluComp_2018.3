#include<stdio.h>
#include<math.h>

#define X 100
#define MR 0.000125

#define L (3.14159265359*2)
#define T 1
#define dX ((double)L)/X

#define EX 2

double f(double u);

double dfdu(double u);

void KT(double u[]);

void Ex1();

void Ex2();

void Ex3();

double minmod(double a, double b);

double sgn(double a);

double minabs(double a, double b);

double max(double a, double b);

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

void KT(double u[]){
  double up[X+3];
  double ux[X+3];
  double uplus[X+3], uminus[X+3];
  double a[X+3];
  double H[X+3];
  double du[X+3];

  double dT = MR*(dX);
  double lbd = MR;
  int t, j;

  //time :: Euler
  //u_j(t+dT) = u_j(t) + ((dU/dt)_j)*dT
  for(t = 0; t*dT <= T; t++){
    //ux vector calculation
    for(j = 1; j <= X+1; j++){
      ux[j] = minmod((u[j]-u[j-1])/(dX), (u[j+1]-u[j])/(dX));
    }

    //uplus and uminus calculation
    for(j = 1; j <= X+1; j++){
      uplus[j] = u[j+1]-((dX)/2)*(ux[j+1]);
      uminus[j] = u[j]+((dX)/2)*(ux[j]);
    }

    //a calculation
    for(j = 1; j <= X+1; j++){
      a[j] = max(dfdu(uplus[j]), dfdu(uminus[j]));
    }

    //H calculation
    for(j = 1; j <= X+1; j++){
      H[j] = (f(uplus[j])+f(uminus[j]))/2 - (a[j]/2)*(uplus[j]-uminus[j]);
    }

    //du calculation
    for(j = 2; j <= X+1; j++){
      du[j] = -(H[j]-H[j-1])/(dX);
    }
    du[1] = du[2];

    //u calculation (Runge-Kutta 2)
    /*for(j = 1; j <= X; j++){
      u[j] = u[j] + ((dT)/2)*(du[j] + dT*du[j+1]);
    }
    u[0] = u[1];
    u[X+2] = u[X+1] = u[X];*/

    //u calculation (Euler)
    for(j = 1; j <= X+1; j++){
      u[j] = u[j] + du[j]*(dT);
    }
    u[0] = u[1];
    u[X+2] = u[X+1];
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

  KT(u);
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

  KT(u);
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

  KT(u);
}

//########### f(u) ###########//
double f(double u){
  switch(EX){
    case 1: return u; break;
    case 2: return pow(u,2)/2; break;
    case 3: return ((pow(u,2)-1)*(pow(u,2)-4))/4; break;
    default: return 0;
  }
}
//###########################//
//########## df/du ##########//
double dfdu(double u){
  switch(EX){
    case 1: return 1; break;
    case 2: return u; break;
    case 3: return pow(u,3)-(5.0/2.0)*u; break;
    default: return 0;
  }
}
//###########################//

//###########################//
//           UTILS           //
//###########################//

double minmod(double a, double b){
  return 0.5*(sgn(a) + sgn(b))*minabs(a,b);
}

double sgn(double a){
  return (0 < a) - (a < 0);
}

double minabs(double a, double b){
  if(fabs(a) < fabs(b)) return fabs(a);
  else return fabs(b);
}

double max(double a, double b){
  if(a > b) return a;
  else return b;
}

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
