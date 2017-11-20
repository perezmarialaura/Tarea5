#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#define Bb 0.2497
#define Bd 5.16
#define Ad 0.3105
#define Ah 64.3 //en kpc
#define N 300 //num datos
#define n 10000 //num iteraciones

double likelihood(double *v_obs, double *v_model);
double *model(double *R_obs, double Mb, double Md, double Mh);
double get_random(void);

int main(void){
  srand(time(NULL)); //sin esto, el rand no cambia al compilar varias veces
  double mb_new, md_new, mh_new, alpha, beta, l_new, l_old;
  float deltab, deltad, deltah;
  double *radius = malloc(N*sizeof(double));
  double *velocities = malloc(N*sizeof(double));
  double *v_old = malloc(N*sizeof(double));
  double *v_new = malloc(N*sizeof(double));
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
    fscanf(data, "%lf %lf \n", &radius[i], &velocities[i]);
  }
  fclose(data);
  //inicializo valores
  mb_walk[0] = get_random()*1000;
  md_walk[0] = get_random()*10000;
  mh_walk[0] = get_random()*100000;
  v_old = model(radius, mb_walk[0], md_walk[0], mh_walk[0]);
  l_walk[0] = likelihood(velocities, v_old);
  deltab = 50;
  deltad = 500;
  deltah = 1000;
  for (i = 0; i < n; i++){
      mb_new = mb_walk[i] + 2*deltab*(get_random()-0.5);
      md_new = md_walk[i] + 2*deltad*(get_random()-0.5);
      mh_new = mh_walk[i] + 2*deltah*(get_random()-0.5);
      v_old = model(radius, mb_walk[i], md_walk[i], mh_walk[i]);
      v_new = model(radius, mb_new, md_new, mh_new);
      l_new = likelihood(velocities, v_new);
      l_old = likelihood(velocities, v_old);
      alpha = l_new/l_old;
      if(alpha >= 1.0){
        mb_walk[i] = mb_new;
        md_walk[i] = md_new;
        mh_walk[i] = mh_new;
        l_walk[i] = l_new;
      }
      else{
        beta = get_random();
        if(alpha >= beta){
          mb_walk[i] = mb_new;
          md_walk[i] = md_new;
          mh_walk[i] = mh_new;
          l_walk[i] = l_new;
        }
        else{
          l_walk[i] = l_old;
        }
      }
    }

  FILE *results =fopen("walks.dat", "w");
  for (i = 0; i < n; i++){
    fprintf(results, "%lf %lf %lf % lf \n", mb_walk[i], md_walk[i], mh_walk[i], l_walk[i]);
  }
    fclose(results);

  return 0;
}

double likelihood(double *v_obs, double *v_model){
  double chi_squared = 0.0;
  int i;
  for(i = 0; i < N; i++){
    chi_squared += 0.5*pow(v_obs[i]-v_model[i],2);
  }
  return exp(-chi_squared);
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
