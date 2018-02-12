#include "particles.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "quadtree.h"
#include "simulations.h"
#include <unistd.h>
#include <sys/time.h>





static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

void simulate(node_t **head, particle_t* particles, particle_t* particles_2, particle_deltas_t* particles_changes, int N, double gravitational_constant, double delta_t, double epsilon, double theta){
  //store particles in struct array
  symplectic_euler(&particles_2, N, gravitational_constant, delta_t, epsilon);
  barnes_hut(head, particles, particles_changes, N, gravitational_constant, delta_t, epsilon, theta);

  # if 1
  //particle_t** particles_ptr, int N, double gravitational_constant, double delta_t, double epsilon

  # else


  # endif
  //update_all(head);



}

void barnes_hut(node_t** head, particle_t* particles_ptr, particle_deltas_t* __dummy, int N, double gravitational_constant, double delta_t, double epsilon, double theta_max){
  //theta_max = 0.4;


  for (int i = 0; i<N; i++){
    double sum_r_x = 0;
    double sum_r_y = 0;

    //traverse the nodes and update the force exertion
    traverse_and_sum(*head, particles_ptr, &sum_r_y, &sum_r_x, theta_max, epsilon, i);
    calculate_deltas(particles_ptr, __dummy, gravitational_constant, sum_r_x, sum_r_y, N, i, delta_t);



  }
  update_all_particles(&particles_ptr, __dummy, N);




}

void calculate_deltas(particle_t* particles, particle_deltas_t* particles_changes, double gravitational_constant,


  double sum_r_x, double sum_r_y, int N, int i, double delta_t){
  //printf("particle pointer1 = %p\n", &particles_changes);
  //subtract own mass?
  //printf("sum b h= %f\n", sum_r_x);
  double acc_i_x_multiplier = -gravitational_constant*sum_r_x*delta_t;
  double acc_i_y_multiplier = -gravitational_constant*sum_r_y*delta_t;
  //
  // particle_t* particle = NULL;
  // particle = particles[i];

  particles_changes[i].x_vel_delta = acc_i_x_multiplier;
  particles_changes[i].y_vel_delta = acc_i_y_multiplier;
  //printf("woow\n");


  double u_n1_x = acc_i_x_multiplier + (particles[i]).x_vel;
  double u_n1_y = acc_i_y_multiplier + (particles[i]).y_vel;
  //printf("or is it down here? \n");
  //store delta position
  particles_changes[i].x_pos_delta = delta_t * u_n1_x;
  particles_changes[i].y_pos_delta = delta_t * u_n1_y;
  //printf("particle pointer 2= %p\n", &particles_changes);
}

void update_all_particles(particle_t** particles, particle_deltas_t* changes, int N){

  for (int i = 0; i<N; i++){
    particle_deltas_t change = changes[i];
    //particle_t particle = particles[i];
    (*particles)[i].x_pos += change.x_pos_delta;
    (*particles)[i].y_pos += change.y_pos_delta;
    (*particles)[i].x_vel += change.x_vel_delta;
    (*particles)[i].y_vel += change.y_vel_delta;
  }
}




void symplectic_euler(particle_t** particles_ptr, int N, double gravitational_constant, double delta_t, double epsilon){
  particle_t *particles = *particles_ptr;

  //struct array storing the particle changes
  struct particle_deltas changes_ptr[N];

  for (int i = 0; i<N; i++){
     particle_t *current_particle;
     current_particle = &particles[i];
     double x_pos = current_particle->x_pos;
     double y_pos = current_particle->y_pos;

     double sum_r_x = 0;
     double sum_r_y = 0;
#if 0
     double r_ij_x[5];
     double r_ij_y[5];

     double r_squared_x[5];
     double r_squared_y[5];

     double r_e[5], r[5];
     double other_mass[5];
     double denominator[5];

    //for loop unrolled
    for(int j=0; j<N; j+=5){


	r_ij_x[0] = (x_pos - other_particle->x_pos);
	r_ij_y[0] = (y_pos - other_particle->y_pos);

	r_ij_x[1] = (x_pos - other_particle->next->x_pos);
        r_ij_y[1] = (y_pos - other_particle->next->y_pos);

	r_ij_x[2] = (x_pos - other_particle->next->next->x_pos);
        r_ij_y[2] = (y_pos - other_particle->next->next->y_pos);

	r_ij_x[3] = (x_pos - other_particle->next->next->next->x_pos);
        r_ij_y[3] = (y_pos - other_particle->next->next->next->y_pos);

	r_ij_x[4] = (x_pos - other_particle->next->next->next->next->x_pos);
        r_ij_y[4] = (y_pos - other_particle->next->next->next->next->y_pos);

	r_squared_x[0] = r_ij_x[0]*r_ij_x[0];
        r_squared_y[0] = r_ij_y[0]*r_ij_y[0];

	r_squared_x[1] = r_ij_x[1]*r_ij_x[1];
        r_squared_y[1] = r_ij_y[1]*r_ij_y[1];

	r_squared_x[2] = r_ij_x[2]*r_ij_x[2];
        r_squared_y[2] = r_ij_y[2]*r_ij_y[2];

	r_squared_x[3] = r_ij_x[3]*r_ij_x[3];
        r_squared_y[3] = r_ij_y[3]*r_ij_y[3];

	r_squared_x[4] = r_ij_x[4]*r_ij_x[4];
        r_squared_y[4] = r_ij_y[4]*r_ij_y[4];

	r[0] = sqrt(r_squared_x[0] + r_squared_y[0]);
	r[1] = sqrt(r_squared_x[1] + r_squared_y[1]);
	r[2] = sqrt(r_squared_x[2] + r_squared_y[2]);
	r[3] = sqrt(r_squared_x[3] + r_squared_y[3]);
	r[4] = sqrt(r_squared_x[4] + r_squared_y[4]);

	r_e[0] = (r[0] + epsilon);
	r_e[1] = (r[1] + epsilon);
	r_e[2] = (r[2] + epsilon);
	r_e[3] = (r[3] + epsilon);
	r_e[4] = (r[4] + epsilon);


        other_mass[0] = other_particle->mass;
	other_mass[1] = other_particle->next->mass;
	other_mass[2] = other_particle->next->next->mass;
	other_mass[3] = other_particle->next->next->next->mass;
	other_mass[4] = other_particle->next->next->next->next->mass;

        //a³ + 3a²b + 3ab² + b³ <- change to this?
        denominator[0] = r_e[0]*r_e[0]*r_e[0];
	denominator[1] = r_e[1]*r_e[1]*r_e[1];
	denominator[2] = r_e[2]*r_e[2]*r_e[2];
	denominator[3] = r_e[3]*r_e[3]*r_e[3];
	denominator[4] = r_e[4]*r_e[4]*r_e[4];
        //double denominator = 7.3;

        sum_r_x += (other_mass[0]*(r_ij_x[0])) / denominator[0];
	sum_r_x += (other_mass[1]*(r_ij_x[1])) / denominator[1];
	sum_r_x += (other_mass[2]*(r_ij_x[2])) / denominator[2];
	sum_r_x += (other_mass[3]*(r_ij_x[3])) / denominator[3];
	sum_r_x += (other_mass[4]*(r_ij_x[4])) / denominator[4];

        sum_r_y += (other_mass[0]*(r_ij_y[0])) / denominator[0];
	sum_r_y += (other_mass[1]*(r_ij_y[1])) / denominator[1];
	sum_r_y += (other_mass[2]*(r_ij_y[2])) / denominator[2];
	sum_r_y += (other_mass[3]*(r_ij_y[3])) / denominator[3];
	sum_r_y += (other_mass[4]*(r_ij_y[4])) / denominator[4];

        //  other_previous_particle = other_particle;
        other_particle = other_particle->next->next->next->next->next;

}
#else
     particle_t other_particle;
     double denominator;
     for (int j = 0; j<N; j++){
          other_particle = particles[j];

          double r_ij_x = (x_pos - other_particle.x_pos);
          double r_ij_y = (y_pos - other_particle.y_pos);
          double r_squared_x = r_ij_x*r_ij_x;
          double r_squared_y = r_ij_y*r_ij_y;
          double r = sqrt(r_squared_x + r_squared_y);
          double r_e = (r + epsilon);
          double other_mass = other_particle.mass;
          denominator = 1/(r_e*r_e*r_e);
          sum_r_x += (other_mass*(r_ij_x)) * denominator;
          sum_r_y += (other_mass*(r_ij_y)) * denominator;

      }
#endif
    //printf("sum simplectic = %f\n", sum_r_x);
     double acc_i_x_multiplier = -gravitational_constant*sum_r_x*delta_t;
     double acc_i_y_multiplier = -gravitational_constant*sum_r_y*delta_t;
    // printf("2 %f \n", sum_r_y);

     changes_ptr[i].x_vel_delta = acc_i_x_multiplier;
     changes_ptr[i].y_vel_delta = acc_i_y_multiplier;

     double u_n1_x = acc_i_x_multiplier + current_particle->x_vel;
     double u_n1_y = acc_i_y_multiplier + current_particle->y_vel;

     //store delta position
     changes_ptr[i].x_pos_delta = delta_t * u_n1_x;
     changes_ptr[i].y_pos_delta = delta_t * u_n1_y;
    }
    //update the new values
    for (int i = 0; i<N; i++){
      particle_deltas_t change = changes_ptr[i];
      //particle_t particle = particles[i];
      (*particles_ptr)[i].x_pos += change.x_pos_delta;
      (*particles_ptr)[i].y_pos += change.y_pos_delta;
      (*particles_ptr)[i].x_vel += change.x_vel_delta;
      (*particles_ptr)[i].y_vel += change.y_vel_delta;
    }
  }



void save(particle_t* particles, int N){
  FILE* file = fopen("results.gal", "wb" );
  for (int i = 0; i<N; i++){
    //printf("%i\n", i);
    particle_t particle = particles[i];
    fwrite(&particle.x_pos, 1, sizeof(double), file );
    fwrite(&particle.y_pos, 1, sizeof(double), file );
    fwrite(&particle.mass, 1, sizeof(double), file );
    fwrite(&particle.x_vel, 1, sizeof(double), file );
    fwrite(&particle.y_vel, 1, sizeof(double), file );
    fwrite(&particle.brightness, 1, sizeof(double), file );
  }
  printf("save complete\n");
}
