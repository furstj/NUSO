#pragma once

#include <cmath>
#include "Grid.hpp"

enum {
  BND_LEFT   = -1,
  BND_RIGHT  = -2,
  BND_TOP    = -3,
  BND_BOTTOM = -4,
};


Grid ChannelGrid(std::size_t Nx, std::size_t Ny); 

