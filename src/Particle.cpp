// #include "Particle.hpp"
// #include "Grid.hpp"
// #include <vector>

// std::vector<particle> initialize_particles(int ppc) {
//   int i, j;
//   std::vector<particle> particles;
//   int num_part = std::sqrt(ppc);
//   for (auto &cell : grid.fluid_cells()) {
//     i = cell->i();
//     j = cell->j();

//     double startx = i * grid.dx() + 0.5 * grid.dx() / num_part;
//     double starty = j * grid.dy() + 0.5 * grid.dy() / num_part;
//     for (int k{0}; k < num_part; k++) {
//       for (int l{0}; l < num_part; l++) {
//         particle p;
//         p.x = startx + k * grid.dx() / num_part;
//         p.y = starty + l * grid.dy() / num_part;

//         particles.push_back(p);
//       }
//     }
//   }
//   return particles;
// }