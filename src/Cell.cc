#include "Cell.h"

Cell::Cell(int level_, double f_, double u_, double v_, double rho_) : 
  level(level_), u(u_), v(v_), rho(rho_)
{
    parent = nullptr;
    for(int i = 0; i < 4; ++i) children[i] = nullptr;
    for(int i = 0; i < 8; ++i) neighbours[i] = nullptr;
    for(int i = 0; i < 9; ++i) f[i] = f_;
}












