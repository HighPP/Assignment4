typedef struct grid {
  double x1, x2, y1, y2;

} grid_t;

typedef struct quad_node {
  double x1, x2, y1, y2;
  double n_particles;
  double total_mass;
  double x_center;
  double y_center;

  struct quad_node* _1;
  struct quad_node* _2;
  struct quad_node* _3;
  struct quad_node* _4;

} node_t;


int create_grid(node_t**, particle_t*, int, int*);
int count_particles(particle_t*, node_t**, int);
void init_tree(node_t**, particle_t*, int);
void print_grid(grid_t*);
int traverse_and_sum(node_t*, particle_t*, double*, double*, double, double, int);
void destroy_tree(node_t* head);
