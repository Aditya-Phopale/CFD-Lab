/*
In this pre-defined file, we define the describe the actual Lid Driven Cavity
grid based on where the fluid lies relative to the cell.
*/
#include "Grid.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#include "Enums.hpp"

Grid::Grid(std::string geom_name, Domain &domain) {
  _domain = domain;

  _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);

  if (geom_name.compare("NONE")) {
    std::vector<std::vector<int>> geometry_data(
        _domain.domain_size_x + 2,
        std::vector<int>(_domain.domain_size_y + 2,
                         0));  //! use domain.size_x instead
    parse_geometry_file(geom_name, geometry_data);
    check_geometry_file(geometry_data);
    assign_cell_types(geometry_data);
    geometry_excluding_ghosts.resize(_domain.size_x,
                                     std::vector<int>(_domain.size_y));
    for (int j = 1; j < _domain.size_y + 1; j++) {
      for (int i = 1; i < _domain.size_x + 1; i++) {
        geometry_excluding_ghosts.at(i - 1).at(j - 1) =
            geometry_data.at(i + domain.imin).at(j + domain.jmin);
      }
    }

  } else {
    build_lid_driven_cavity();
  }
}

void Grid::build_lid_driven_cavity() {
  std::vector<std::vector<int>> geometry_data(
      _domain.domain_size_x + 2,
      std::vector<int>(_domain.domain_size_y + 2, 0));

  for (int i = _domain.imin; i < _domain.imax; ++i) {
    for (int j = _domain.jmin; j < _domain.jmax; ++j) {
      // Bottom, left and right walls: no-slip
      if (i == 0 || j == 0 || i == _domain.domain_size_x + 1) {
        geometry_data.at(i).at(j) = LidDrivenCavity::fixed_wall_id;
      }
      // Top wall: moving wall
      else if (j == _domain.domain_size_y + 1) {
        geometry_data.at(i).at(j) = LidDrivenCavity::moving_wall_id;
      }
    }
  }

  assign_cell_types(geometry_data);

  geometry_excluding_ghosts.resize(_domain.domain_size_x,
                                   std::vector<int>(_domain.domain_size_y));
  for (int j = _domain.jmin; j < jmax(); j++) {
    for (int i = _domain.imin; i < imax(); i++) {
      geometry_excluding_ghosts.at(i).at(j) = geometry_data.at(i + 1).at(j + 1);
    }
  }
}

void Grid::assign_cell_types(std::vector<std::vector<int>> &geometry_data) {
  int i = 0;
  int j = 0;
  for (int j_geom = _domain.jmin; j_geom < _domain.jmax; ++j_geom) {
    { i = 0; }
    for (int i_geom = _domain.imin; i_geom < _domain.imax; ++i_geom) {
      if (geometry_data.at(i_geom).at(j_geom) == cellID::fluid) {
        _cells(i, j) = Cell(i, j, cell_type::FLUID);
        if (i_geom == _domain.imin || i_geom == _domain.imax - 1 ||
            j_geom == _domain.jmin || j_geom == _domain.jmax - 1) {
          _buffer.push_back(&_cells(i, j));
        } else {
          _fluid_cells.push_back(&_cells(i, j));
        }

      } else if (geometry_data.at(i_geom).at(j_geom) == cellID::fixed_wall_3) {
        _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL3,
                            geometry_data.at(i_geom).at(j_geom));

        _fixed_wall_cells.push_back(&_cells(i, j));

      } else if (geometry_data.at(i_geom).at(j_geom) == cellID::fixed_wall_4) {
        _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL4,
                            geometry_data.at(i_geom).at(j_geom));

        _fixed_wall_cells.push_back(&_cells(i, j));

      } else if (geometry_data.at(i_geom).at(j_geom) == cellID::fixed_wall_5) {
        _cells(i, j) = Cell(i, j, cell_type::ADIABATIC_WALL,
                            geometry_data.at(i_geom).at(j_geom));
        _adiabatic_cells.push_back(&_cells(i, j));

      } else if (geometry_data.at(i_geom).at(j_geom) == cellID::inflow) {
        _cells(i, j) =
            Cell(i, j, cell_type::INLET, geometry_data.at(i_geom).at(j_geom));
        _inlet_cells.push_back(&_cells(i, j));

      } else if (geometry_data.at(i_geom).at(j_geom) == cellID::outflow) {
        _cells(i, j) =
            Cell(i, j, cell_type::OUTLET, geometry_data.at(i_geom).at(j_geom));
        _outlet_cells.push_back(&_cells(i, j));

      } else if (geometry_data.at(i_geom).at(j_geom) ==
                 LidDrivenCavity::moving_wall_id) {
        _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL,
                            geometry_data.at(i_geom).at(j_geom));

        _moving_wall_cells.push_back(&_cells(i, j));
      } else if (geometry_data.at(i_geom).at(j_geom) == cellID::empty) {
        _cells(i, j) =
            Cell(i, j, cell_type::EMPTY, geometry_data.at(i_geom).at(j_geom));
      }
      ++i;
    }
    ++j;
  }

  // Corner cell neighbour assigment
  // Bottom-Left Corner
  i = 0;
  j = 0;
  _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
  _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
  if (_cells(i, j).neighbour(border_position::TOP)->type() ==
      cell_type::FLUID) {
    _cells(i, j).add_border(border_position::TOP);
  }
  if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
      cell_type::FLUID) {
    _cells(i, j).add_border(border_position::RIGHT);
  }
  // Top-Left Corner
  i = 0;
  j = _domain.size_y + 1;
  _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
  _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
  if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
      cell_type::FLUID) {
    _cells(i, j).add_border(border_position::BOTTOM);
  }
  if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
      cell_type::FLUID) {
    _cells(i, j).add_border(border_position::RIGHT);
  }

  // // Top-Right Corner
  i = _domain.size_x + 1;
  j = Grid::_domain.size_y + 1;
  _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
  _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
  if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
      cell_type::FLUID) {
    _cells(i, j).add_border(border_position::BOTTOM);
  }
  if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
      cell_type::FLUID) {
    _cells(i, j).add_border(border_position::LEFT);
  }

  // Bottom-Right Corner
  i = Grid::_domain.size_x + 1;
  j = 0;
  _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
  _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
  if (_cells(i, j).neighbour(border_position::TOP)->type() ==
      cell_type::FLUID) {
    _cells(i, j).add_border(border_position::TOP);
  }
  if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
      cell_type::FLUID) {
    _cells(i, j).add_border(border_position::LEFT);
  }
  // Bottom cells
  j = 0;
  for (int i = 1; i < _domain.size_x + 1; ++i) {
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    if (_cells(i, j).type() != cell_type::FLUID) {
      if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::RIGHT);
      }
      if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::LEFT);
      }
      if (_cells(i, j).neighbour(border_position::TOP)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::TOP);
      }
    }
  }

  // Top Cells
  j = Grid::_domain.size_y + 1;

  for (int i = 1; i < _domain.size_x + 1; ++i) {
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    if (_cells(i, j).type() != cell_type::FLUID) {
      if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::RIGHT);
      }
      if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::LEFT);
      }
      if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::BOTTOM);
      }
    }
  }

  // Left Cells
  i = 0;
  for (int j = 1; j < _domain.size_y + 1; ++j) {
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    if (_cells(i, j).type() != cell_type::FLUID) {
      if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::RIGHT);
      }
      if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::BOTTOM);
      }
      if (_cells(i, j).neighbour(border_position::TOP)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::TOP);
      }
    }
  }
  // Right Cells
  i = Grid::_domain.size_x + 1;
  for (int j = 1; j < _domain.size_y + 1; ++j) {
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    if (_cells(i, j).type() != cell_type::FLUID) {
      if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::LEFT);
      }
      if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::BOTTOM);
      }
      if (_cells(i, j).neighbour(border_position::TOP)->type() ==
          cell_type::FLUID) {
        _cells(i, j).add_border(border_position::TOP);
      }
    }
  }

  // Inner cells
  for (int i = 1; i < _domain.size_x + 1; ++i) {
    for (int j = 1; j < _domain.size_y + 1; ++j) {
      _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
      _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
      _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
      _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);

      if (_cells(i, j).type() != cell_type::FLUID &&
          _cells(i, j).type() != cell_type::EMPTY) {
        if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::FLUID &&
            _cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::FLUID) {
          _cells(i, j).add_border(border_position::NORTHWEST);
          continue;
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::FLUID &&
            _cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::FLUID) {
          _cells(i, j).add_border(border_position::NORTHEAST);
          continue;
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::FLUID &&
            _cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::FLUID) {
          _cells(i, j).add_border(border_position::SOUTHEAST);
          continue;
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::FLUID &&
            _cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::FLUID) {
          _cells(i, j).add_border(border_position::SOUTHWEST);
          continue;
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
            cell_type::FLUID) {
          _cells(i, j).add_border(border_position::LEFT);
          continue;
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
            cell_type::FLUID) {
          _cells(i, j).add_border(border_position::RIGHT);
          continue;
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
            cell_type::FLUID) {
          _cells(i, j).add_border(border_position::BOTTOM);
          continue;
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() ==
            cell_type::FLUID) {
          _cells(i, j).add_border(border_position::TOP);
          continue;
        }
      }

      if (_cells(i, j).type() == cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::NORTHEASTWESTSOUTH);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::EASTWESTSOUTH);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::NORTHWESTSOUTH);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::NORTHEASTSOUTH);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::NORTHEASTWEST);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::EASTWEST);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::NORTHSOUTH);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::SOUTHWEST);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::NORTHEAST);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::BOTTOM)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::SOUTHEAST);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() ==
                cell_type::EMPTY &&
            _cells(i, j).neighbour(border_position::LEFT)->type() ==
                cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::NORTHWEST);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() ==
            cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::RIGHT);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() ==
            cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::LEFT);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() ==
            cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::TOP);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() ==
            cell_type::EMPTY) {
          _cells(i, j).add_border(border_position::BOTTOM);
          _surface_cells.push_back(&_cells(i, j));
          continue;
        }
      }
    }
  }
}

void Grid::parse_geometry_file(std::string filedoc,
                               std::vector<std::vector<int>> &geometry_data) {
  int numcols, numrows, depth;

  std::ifstream infile(filedoc);
  std::stringstream ss;
  std::string inputLine = "";

  // First line : version
  getline(infile, inputLine);
  if (inputLine.compare("P2") != 0) {
    std::cerr << "First line of the PGM file should be P2" << std::endl;
  }

  // Second line : comment
  getline(infile, inputLine);

  // Continue with a stringstream
  ss << infile.rdbuf();
  // Third line : size
  ss >> numrows >> numcols;
  // Fourth line : depth
  ss >> depth;

  // Following lines : data
  for (int col = numcols - 1; col > -1; --col) {
    for (int row = 0; row < numrows; ++row) {
      ss >> geometry_data[row][col];
    }
  }

  infile.close();
}

void Grid::check_geometry_file(std::vector<std::vector<int>> &geometry_data) {
  for (int j = _domain.jmin + 1; j < _domain.jmax - 1; ++j) {
    for (int i = _domain.imin + 1; i < _domain.imax - 1; ++i) {
      if (geometry_data.at(i).at(j) == 3 || geometry_data.at(i).at(j) == 4 ||
          geometry_data.at(i).at(j) == 5) {
        int sum = geometry_data.at(i + 1).at(j) +
                  geometry_data.at(i - 1).at(j) +
                  geometry_data.at(i).at(j + 1) + geometry_data.at(i).at(j - 1);

        if (sum <= 5) {
          geometry_data.at(i).at(j) = 0;
        }
      }
    }
  }
}

int Grid::imax() const { return _domain.size_x; }
int Grid::jmax() const { return _domain.size_y; }

int Grid::imaxb() const { return _domain.size_x + 2; }
int Grid::jmaxb() const { return _domain.size_y + 2; }

Cell Grid::cell(int i, int j) const { return _cells(i, j); }

double Grid::dx() const { return _domain.dx; }

double Grid::dy() const { return _domain.dy; }

const Domain &Grid::domain() const { return _domain; }

std::vector<Cell *> &Grid::fluid_cells() { return _fluid_cells; }

void Grid::set_particles(int ppc) {
  // _particles = initialize_particles(ppc);
  // for (auto &elem : p) {
  //   _particles.push_back(elem);
  // }

  int i, j;
  //std::vector<particle> particles;
  int num_part = std::sqrt(ppc);
  for (auto &cell : _fluid_cells) {
    i = cell->i();
    j = cell->j();

    double startx = i * _dx + 0.5 * _dx / num_part;
    double starty = j * _dy + 0.5 * _dy / num_part;
    for (int k{0}; k < num_part; k++) {
      for (int l{0}; l < num_part; l++) {
        particle p;
        p.x = startx + k * _dx / num_part;
        p.y = starty + l * _dy / num_part;

        _particles.push_back(p);
      }
    }
  }
  //return particles;
}

const std::vector<Cell *> &Grid::fixed_wall_cells() const {
  return _fixed_wall_cells;
}

const std::vector<Cell *> &Grid::moving_wall_cells() const {
  return _moving_wall_cells;
}

const std::vector<Cell *> &Grid::inlet_cells() const { return _inlet_cells; }

const std::vector<Cell *> &Grid::outlet_cells() const { return _outlet_cells; }

const std::vector<Cell *> &Grid::adiabatic_cells() const {
  return _adiabatic_cells;
}

std::vector<Cell *> &Grid::surface_cells(){
  return _surface_cells;
}

const std::vector<Cell *> &Grid::buffer() const { return _buffer; }

const std::vector<std::vector<int>> &Grid::get_geometry_excluding_ghosts()
    const {
  return geometry_excluding_ghosts;
}
