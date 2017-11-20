#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#define Bb 0.2497
#define Bd 5.16
#define Ad 0.3105
#define Ah 64.3 //en kpc
#define N 300 //num datos
#define n 100000 //num iteraciones

double likelihood(double *v_obs, double *v_model);
double *model(double *R_obs, double Mb, double Md, double Mh);
double get_random(void);

int main(void){
  srand(time(NULL)); //sin esto, el rand no cambia al compilar varias veces
  double l_best, mb_best, md_best, mh_best, mb_prime, md_prime, mh_prime, alpha, beta, l_prime, l_init, deltab, deltad, deltah;
  double *r_obs = malloc(N*sizeof(double));
  double *v_obs = malloc(N*sizeof(double));
  double *v_init = malloc(N*sizeof(double));
  double *v_prime = malloc(N*sizeof(double));
  double *mb_walk = malloc(n*sizeof(double));
  double *md_walk = malloc(n*sizeof(double));
  double *mh_walk = malloc(n*sizeof(double));
  double *l_walk = malloc(n*sizeof(double));

  FILE *data = fopen("RadialVelocities.dat", "r");
  int i;
  //elimino la primera fila
  char buffer[100];
  fgets(buffer, 100, data);
  //e importo los datos
  for (i = 0; i < N; i++){
    fscanf(data, "%lf %lf \n", &r_obs[i], &v_obs[i]);
  }
  fclose(data);
  //inicializo valores
  mb_walk[0] = 300; //Para comenzar cercano, se estiman los iniciales como un primer guess
  md_walk[0] = 7000;
  mh_walk[0] = 10000;

  v_init = model(r_obs, mb_walk[0], md_walk[0], mh_walk[0]);
  l_walk[0] = likelihood(v_obs, v_init);
  deltab = 0.05; //como tienen diferentes ordenes deberían variar los delta
  deltad = 1;
  deltah = 5;
  for (i = 1; i < n; i++){ //Algoritmo de MonteCarlo
      mb_prime = mb_walk[i-1] + 2*deltab*(get_random()-0.5);
      md_prime = md_walk[i-1] + 2*deltad*(get_random()-0.5);
      mh_prime = mh_walk[i-1] + 2*deltah*(get_random()-0.5);

      v_init = model(r_obs, mb_walk[i-1], md_walk[i-1], mh_walk[i-1]);
      v_prime = model(r_obs, mb_prime, md_prime, mh_prime);

      l_prime = likelihood(v_obs, v_prime);
      l_init = likelihood(v_obs, v_init);

      alpha = l_prime/l_init;
      if(alpha >= 1.0){
        mb_walk[i] = mb_prime;
        md_walk[i] = md_prime;
        mh_walk[i] = mh_prime;
        l_walk[i] = l_prime;
      }
      else{
        beta = get_random();
        if(alpha >= beta){
          mb_walk[i] = mb_prime;
          md_walk[i] = md_prime;
          mh_walk[i] = mh_prime;
          l_walk[i] = l_prime;
        }
        else{
          mb_walk[i] = mb_walk[i-1];
          md_walk[i] = md_walk[i-1];
          mh_walk[i] = mh_walk[i-1];
          l_walk[i] = l_init;
        }
      }
    }
    //busco el mayor likelihood y sus M correspondientes
l_best = l_walk[0];
for (i = 1; i < n; i++){
  if (l_walk[i] > l_best){
    l_best = l_walk[i];
    mb_best = mb_walk[i];
    md_best = md_walk[i];
    mh_best = mh_walk[i];
  }
}
//impresión en consola
printf("%s %f %s %f %s %f \n", "Mb= ",mb_best, ", Md= ", md_best, ", Mh = ", mh_best);
//.dat que usará python
FILE *results = fopen("M.dat","w");
fprintf(results, "%f %f %f \n", mb_best, md_best, mh_best);
fclose(results);

  return 0;
}

double likelihood(double *v_obs, double *v_model){
  double chi_squared = 0.0;
  int i;
  for(i = 0; i < N; i++){
    chi_squared += 0.5*pow(v_obs[i]-v_model[i],2);
  }
  return exp(-chi_squared/1000.);
}

double *model(double *R_obs, double Mb, double Md, double Mh){
  int i;
  double *V = malloc(N*sizeof(double));
  for (i = 0; i < N; i++){
    V[i] = (R_obs[i]*pow(Mb, 0.5)/pow(pow(R_obs[i],2)+pow(Bb,2),0.75)) + (R_obs[i]*pow(Mb,0.5)/pow(pow(R_obs[i],2)+pow(Bd+Ad,2),0.75)) + (pow(Mh,0.5)/pow(pow(R_obs[i],2) + pow(Ah,2),0.25));
  }
  return V;
}

double get_random(void){
  return (double) rand()/RAND_MAX;
}
