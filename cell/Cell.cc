#include "Cell.h"

Cell::Cell(int level_, double u_, double v_, double rho_) : 
  level(level_), u(u_), v(v_), rho(rho_)
{
    parent = nullptr;
    for(int i = 0; i < 4; ++i) children[i] = nullptr;
    e = nullptr;
    ne = nullptr;
    n = nullptr;
    nw = nullptr;
    w = nullptr;
    sw = nullptr;
    s = nullptr;
    se = nullptr;
    f = 0;
    fe = 0;
    fne = 0;
    fn = 0;
    fnw = 0;
    fw = 0;
    fsw = 0;
    fs = 0;
    fse = 0;
    fbuffer = 0;
    // for(int i = 0; i < 8; ++i) neighbours[i] = nullptr;
    // for(int i = 0; i < 9; ++i) f[i] = 0;
    // Iniitalize f via equilibrium distribution function.
}

double Cell::get_velocity_magnitude()
{
  return sqrt(u*u + v*v);
}










