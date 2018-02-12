typedef struct particle {
  double x_pos, y_pos, mass, x_vel, y_vel, brightness;

} particle_t;

typedef struct particle_deltas {
  double x_pos_delta, y_pos_delta, x_vel_delta, y_vel_delta;

} particle_deltas_t;
