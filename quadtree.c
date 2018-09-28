#include "particles.h"
#include "quadtree.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <assert.h>
#include <string.h>
#include <math.h>

int create_grid(node_t** root, particle_t* particles, int N, int* leaves){
  double x1, x2, y1, y2;
  x1 = (*root)->x1;
  x2 = (*root)->x2;
  y1 = (*root)->y1;
  y2 = (*root)->y2;

  int counter = 0;
  node_t** quads = NULL;
  quads = malloc(4 * sizeof(node_t*));

  double y_lower = (y1 + y2)*0.5;
  double y_upper = y2;

  //Maybe we can cobine together the two for cycles unrolling the second one like this:
  
  /*
    for (int i = 0; i <2; i++){
    double x_lower = x1;
    double x_upper = (x2 + x1)*0.5;
    
    quads[0+i] = malloc(sizeof(node_t));
    quads[0+i]->x1 = x_lower;
    quads[0+i]->x2 = x_upper;
    quads[0+i]->y1 = y_lower;
    quads[0+i]->y2 = y_upper;
    quads[0+i]->x_center = 0;
    quads[0+i]->y_center = 0;
    quads[0+i]->total_mass = 0;
    quads[0+i]->n_particles = 0;
    quads[0+i]->_1 = NULL;
    quads[0+i]->_2 = NULL;
    quads[0+i]->_3 = NULL;
    quads[0+i]->_4 = NULL;

    x_lower = x_upper;
    x_upper = x2;
      
    quads[1+i] = malloc(sizeof(node_t));
    quads[1+i]->x1 = x_lower;
    quads[1+i]->x2 = x_upper;
    quads[1+i]->y1 = y_lower;
    quads[1+i]->y2 = y_upper;
    quads[1+i]->x_center = 0;
    quads[1+i]->y_center = 0;
    quads[1+i]->total_mass = 0;
    quads[1+i]->n_particles = 0;
    quads[1+i]->_1 = NULL;
    quads[1+i]->_2 = NULL;
    quads[1+i]->_3 = NULL;
    quads[1+i]->_4 = NULL;
      
    //adjust y-traverse
    y_upper = y_lower;
    y_lower = y1;
  }
  */
  for (int i = 0; i <2; i++){
    double x_lower = x1;
    double x_upper = (x2 + x1)*0.5;
    for (int j = 0; j<2; j++){
      //Instead of using counter should we just use the variable j?
      quads[counter] = malloc(sizeof(node_t));
      quads[counter]->x1 = x_lower;
      quads[counter]->x2 = x_upper;
      quads[counter]->y1 = y_lower;
      quads[counter]->y2 = y_upper;
      quads[counter]->x_center = 0;
      quads[counter]->y_center = 0;
      quads[counter]->total_mass = 0;
      quads[counter]->n_particles = 0;
      quads[counter]->_1 = NULL;
      quads[counter]->_2 = NULL;
      quads[counter]->_3 = NULL;
      quads[counter]->_4 = NULL;

      x_lower = x_upper;
      x_upper = x2;
      counter ++;
    }
    //adjust y-traverse
    y_upper = y_lower;
    y_lower = y1;
  }
  
  //Here we already divide the main quad in 4, should we count the particle here and if is bigger than one divide tha quad otherwise 
  //point to NULL? 
  (*root)->_1 = quads[0];
  (*root)->_2 = quads[1];
  (*root)->_3 = quads[2];
  (*root)->_4 = quads[3];

  //It can return an array of particles that we're gonna check for the next children
  int no_particles = count_particles(particles, root, N);
  
  //Should we end the function here and return the n. of particles and recall the function until the number of particles is = 0?
  if (no_particles == 1){
    //only for assertion purposes - make sure all levaes been created
    (*leaves)++;
  }
 
  //If a leaf contains more than one particle then split it again and check how many particles are now inside
  //Are we splitting again all 4? We should split only the leaf that contains more than 1 particle while leaving the other intact
  if (no_particles > 1){
    int n_child_particles = 0;
    
    //Should we check how many particles does the function return? 
    //Here are we dividing again in 4 each quad but we have checked the number of particle in the root? Maybe not 
    create_grid(&quads[0], particles, N, leaves);
    create_grid(&quads[1], particles, N, leaves);
    create_grid(&quads[2], particles, N, leaves);
    create_grid(&quads[3], particles, N, leaves);
    n_child_particles += quads[0]->n_particles;
    n_child_particles += quads[1]->n_particles;
    n_child_particles += quads[2]->n_particles;
    n_child_particles += quads[3]->n_particles;
    assert(no_particles == n_child_particles);
  }
  free(quads);
  return no_particles;
}

int count_particles(particle_t* particles, node_t** node, int N){
  double x1, x2, y1, y2;
  x1 = (*node)->x1;
  x2 = (*node)->x2;
  y1 = (*node)->y1;
  y2 = (*node)->y2;
  int counter = 0;
  double x_sum = 0;
  double y_sum = 0;
  double mass = 0;

  //optimzie with register comparisons?
  //I thought about changing the struct that we are pointing at in each step so each time we don't have to go through all the particles
  //again but only to the ones contained in the previous bigger quad, then the research should go faster. (1)
  for (int i = 0; i<N; i++){
    if (particles[i].x_pos >= x1 && particles[i].x_pos <= x2 && particles[i].y_pos > y1 && particles[i].y_pos <= y2){
      //save these particles in a new particle_t structure and pass this to the child node so we're gonna check for the amount of particles
      //through a smaller array of structs 
      counter ++;
      x_sum += particles[i].x_pos;
      y_sum += particles[i].y_pos;
      mass += particles[i].mass;
    }
  }
  (*node)->total_mass = mass;
  (*node)->x_center = x_sum/counter;
  (*node)->y_center = y_sum/counter;
  (*node)->n_particles = counter;

  return counter;
}

/*
particle_t count_particles(particle_t* particles, node_t** node, int N){
  particle_t* in_particles;
  double x1, x2, y1, y2;
  x1 = (*node)->x1;
  x2 = (*node)->x2;
  y1 = (*node)->y1;
  y2 = (*node)->y2;
  int counter = 0;
  double x_sum = 0;
  double y_sum = 0;
  double mass = 0;

  in_particles =(particle_t*)malloc(sizeof(particle_t));
  //optimzie with register comparisons?
  //I thought about changing the struct that we are pointing at in each step so each time we don't have to go through all the particles
  //again but only to the ones contained in the previous bigger quad, then the research should go faster. (1)
  for (int i = 0; i<N; i++){
    if (particles[i].x_pos >= x1 && particles[i].x_pos <= x2 && particles[i].y_pos > y1 && particles[i].y_pos <= y2){
      //save these particles in a new particle_t structure and pass this to the child node so we're gonna check for the amount of particles
      //through a smaller array of structs 
      counter ++;
      x_sum += particles[i].x_pos;
      y_sum += particles[i].y_pos;
      mass += particles[i].mass;
      
      realloc(in_particle, count*sizeof(particle_t))
      memcpy(in_particles[count], particles[i], sizeof(particle_t)); //Does it make sense?
    }
  }
  (*node)->total_mass = mass;
  (*node)->x_center = x_sum/counter;
  (*node)->y_center = y_sum/counter;
  (*node)->n_particles = counter;

  return in_particles;
}
*/
void init_tree(node_t**head, particle_t* particles, int N){
  (*head)->x1 = 0; (*head)->x2 = 1; (*head)->y1 = 0; (*head)->y2 = 1;
  int n_leaves = 0;
  create_grid(head, particles, N, &n_leaves);
}


int traverse_and_sum(node_t* quadrant, particle_t* particle_ptr, double* sum_r_y,
  double* sum_r_x, double theta_max, double epsilon, int i){

  //printf("traversing\n");

  double width = quadrant->x2 - quadrant->x1;
  double r_ij_x = (particle_ptr[i]).x_pos - quadrant->x_center;
  // printf("hej b h%f \n", r_ij_x);
  double r_ij_y = (particle_ptr[i]).y_pos - quadrant->y_center;

  double r_squared_x = r_ij_x*r_ij_x;
  double r_squared_y = r_ij_y*r_ij_y;
  double r = sqrt(r_squared_x + r_squared_y);

  double theta = width/r;
  // printf("theta = %f particles = %f x center = %f\n", theta, quadrant->n_particles, quadrant->total_mass);
  // printf("theta_max = %f\n", theta_max);

  if (theta<theta_max || quadrant->n_particles==1){
    //printf("n = %f\n", quadrant->n_particles);
    double r_e = (r + epsilon);
    double other_mass = quadrant->total_mass;
    double denominator = 1/(r_e*r_e*r_e);
    *sum_r_x += (other_mass*(r_ij_x)) * denominator;
    *sum_r_y += (other_mass*(r_ij_y)) * denominator;
    return 1;
  }
    if (quadrant->_1){
      //printf("quadrant 1\n");
      traverse_and_sum(quadrant->_1, particle_ptr, sum_r_y, sum_r_x, theta_max, epsilon, i);
    } if(quadrant->_2){
    //  printf("quadrant 2\n");
      traverse_and_sum(quadrant->_2, particle_ptr, sum_r_y, sum_r_x, theta_max, epsilon, i);
    } if (quadrant->_3){
    //  printf("quadrant 3\n");
      traverse_and_sum(quadrant->_3, particle_ptr, sum_r_y, sum_r_x, theta_max, epsilon, i);
    } if (quadrant->_4){
    //  printf("quadrant 3\n");
      traverse_and_sum(quadrant->_4, particle_ptr, sum_r_y, sum_r_x, theta_max, epsilon, i);
    }
    return 1;
}

void destroy_tree(node_t* quadrant){
  if (quadrant->_1){
    destroy_tree(quadrant->_1);
  } if(quadrant->_2){
    destroy_tree(quadrant->_2);
  } if (quadrant->_3){
    destroy_tree(quadrant->_3);
  } if (quadrant->_4){
    destroy_tree(quadrant->_4);
  }
  free(quadrant);
}
