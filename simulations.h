void simulate(node_t**, particle_t*, particle_t*, particle_deltas_t*, int, double, double, double, double);
void update_all_particles(particle_t**, particle_deltas_t*, int);


void symplectic_euler(particle_t**, int, double, double, double);
void barnes_hut(node_t**, particle_t*, particle_deltas_t*, int, double, double delta_t, double, double);
void calculate_deltas(particle_t*, particle_deltas_t*, double, double , double, int, int, double);

void save(particle_t*, int);
