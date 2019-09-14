#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "particles.h"
#include "quadtree.h"
#include "simulations.h"
#include <unistd.h>
#include <sys/time.h>
#include <assert.h>
//
// void simulate();
// void save();

static double get_wall_seconds();

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

int main(int argc, char **argv){

  double N, n_steps; //graphics;
  char filename[255];
  double delta_t;
  double theta_max;

  if (argc != 7){
    printf("Please enter 6 arguments (N file name nsteps delta_t theta_max graphics)");
    return -1;
  } else {
    //params
    N = atoi(argv[1]);
    strcpy(filename, argv[2]);
    n_steps = atoi(argv[3]);
    delta_t = atof(argv[4]);
    theta_max = atof(argv[5]);
    //graphics = atoi(argv[5]);

  }
///This is just a try to see how Github works 
  double gravitational_constant = 100/N;
  printf("grav const = %f\n", gravitational_constant);
  double epsilon = 1e-3;
  delta_t = 1e-5;

  printf("dt = %f\n", delta_t);
  printf("eps = %f\n", epsilon);
  // n_steps = 100;

  double time1;
  double time2;

  double x_pos;
  double y_pos;
  double mass;
  double x_vel;
  double y_vel;
  double dummy;

  //Stor particles on stack
  struct particle particles[(int) N];
  struct particle particles_2[(int) N];
  struct particle_deltas particles_changes[(int) N];

  //read data from byte file
  time1 = get_wall_seconds();
  FILE* fp1;
  fp1 = fopen(filename, "rb");
  char *buffer;
  size_t bytes_read;
  for (int i = 0; i<N; i++){
    buffer = (char *)malloc(sizeof(double));

    bytes_read = fread(buffer, sizeof(double), 1, fp1);
    if (bytes_read == 0){
      printf("You've reached the end of the file, stopping here...\n");
      sleep(2);
      return -1;
    }


    memcpy(&x_pos, buffer,sizeof(double));
    memset(buffer, 0, sizeof(double));
    particles[i].x_pos = x_pos;
    particles_2[i].x_pos = x_pos;

    bytes_read = fread(buffer, sizeof(double), 1, fp1);
    memcpy(&y_pos, buffer,sizeof(double));
    memset(buffer, 0, sizeof(double));
    particles[i].y_pos = y_pos;
    particles_2[i].y_pos = y_pos;

    bytes_read = fread(buffer, sizeof(double), 1, fp1);
    memcpy(&mass, buffer,sizeof(double));
    memset(buffer, 0, sizeof(double));
    particles[i].mass = mass;
    particles_2[i].mass = mass;

    bytes_read = fread(buffer, sizeof(double), 1, fp1);
    memcpy(&x_vel, buffer,sizeof(double));
    memset(buffer, 0, sizeof(double));
    particles[i].x_vel = x_vel;
    particles_2[i].x_vel = x_vel;


    bytes_read = fread(buffer, sizeof(double), 1, fp1);
    memcpy(&y_vel, buffer,sizeof(double));
    memset(buffer, 0, sizeof(double));
    particles[i].y_vel = y_vel;
    particles_2[i].y_vel = y_vel;

    bytes_read = fread(buffer, sizeof(double), 1, fp1);
    memcpy(&dummy, buffer,sizeof(double));
    memset(buffer, 0, sizeof(double));
    particles[i].brightness = dummy;
    particles_2[i].brightness = dummy;
    //printf("particle data = %f\n", current_particle->y_vel);

    particles_changes[i].x_pos_delta = 0.0;
    particles_changes[i].y_pos_delta = 0.0;
    particles_changes[i].x_vel_delta = 0.0;
    particles_changes[i].y_vel_delta = 0.0;
    // printf("particles data x/y = [%f %f]\n", x_pos, y_pos);
    free(buffer);

  }
  fclose(fp1);
  printf("reading took %f seconds  \n", get_wall_seconds()-time1);
;



  time2 = get_wall_seconds();
  int step = 0;
  //double theta = 0.2;
  while (step < n_steps){

    //maybe head doesn't have to be redfined
    node_t* head = NULL;
    head = malloc(sizeof(node_t));
    init_tree(&head, particles, N);
    simulate(&head, particles, particles_2, particles_changes, N, gravitational_constant, delta_t, epsilon, theta_max);
    printf("with barnes_hut %f; ", particles[1].x_pos);
    printf("with symplectic %f; \n", particles_2[1].x_pos);
    destroy_tree(head);
    step++;
  }

  printf("simulating took %f seconds \n", get_wall_seconds()-time2);
  save(particles, N);
  //print_linked_list(head);
  return 0;
}
