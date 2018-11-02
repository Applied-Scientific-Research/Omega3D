#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <variant>
#define _USE_MATH_DEFINES
#include <cmath>

const size_t Dimensions = 3;

// templated solution type

// templated on the accumulator precision
//template <class A>
//class Solver {
//public:
//private:
//};

enum solution_t {
  direct_cpu   = 1,
  direct_vc    = 2,
  direct_glsl  = 3,
  treecode_cpu = 4
};

enum elem_t {
  active   = 1,
  reactive = 2,
  inert    = 3
};

enum move_t {
  lagrangian = 1,
  bodybound  = 2,
  fixed      = 3
};


// the superclass


template <class S>
class ElementBase {
public:
  ElementBase<S>(const size_t _n, const elem_t _e, const move_t _m) :
      n(_n), E(_e), M(_m) {
  }

  size_t getn() const { return n; }
  const std::vector<S>& get_pos() const { return x; }
  const std::vector<S>& get_rad() const { return r; }
  const std::vector<S>& get_str() const { return *s; }
  std::vector<S>& get_vel() { return u; }

  void zero_vels() { for( S& val : u ) val = 0.0; }
  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    for (size_t i=0; i<getn(); ++i) {
      u[3*i+0] = _fs[0] + u[3*i+0] * 0.25/M_PI;
      u[3*i+1] = _fs[1] + u[3*i+1] * 0.25/M_PI;
      u[3*i+2] = _fs[2] + u[3*i+2] * 0.25/M_PI;
    }
  }
  void move(const double _dt) {
    if (M == lagrangian) {
      std::cout << "  Moving" << to_string() << std::endl;

      // update positions
      for (size_t i=0; i<3*n; ++i) {
        x[i] += (S)_dt * u[i];
      }

      // update elongation
      // update strengths
    }
  }

  std::string to_string() const {
    std::string mystr = " " + std::to_string(n);
    if (E == active) {
      mystr += " Active";
    } else if (E == reactive) {
      mystr += " Reactive";
    } else {
      mystr += " Inert";
    }
    if (M == lagrangian) {
      mystr += " Lagrangian";
    } else if (M == bodybound) {
      mystr += " Body-fixed";
    } else {
      mystr += " Fixed";
    }
    return mystr;
  }

protected:
  // active, reactive, or inert?
  elem_t E;
  // how does it move? use move_t or Body*
  move_t M;
  //Body* b = nullptr;
  //Move_t get_move_type() {
  //}

  // common arrays for all derived types
  size_t n;
  // state vector
  std::vector<S> x;            // position
  std::vector<S> r;            // thickness/radius
  std::optional<std::vector<S>> s;    // strength
  // time derivative of state vector
  std::vector<S> u;            // velocity
  std::optional<std::vector<S>> dsdt;    // strength change
};



// two types of collections (for now)

// 0-D elements
template <class S>
class Points: public ElementBase<S> {
public:
  Points(const size_t _n, const elem_t _e, const move_t _m) :
      ElementBase<S>(_n, _e, _m) {
    // this initialization specific to Points
    this->x.resize(3*_n);
    this->r.resize(_n);
    this->elong.resize(_n);
    for (size_t i=0; i<_n; ++i) {
      this->x[3*i+0] = -1.0 + 2.0*(S)rand()/(S)RAND_MAX;
      this->x[3*i+1] = -1.0 + 2.0*(S)rand()/(S)RAND_MAX;
      this->x[3*i+2] = -1.0 + 2.0*(S)rand()/(S)RAND_MAX;
      this->r[i] = 0.01;
      this->elong[i] = 1.0;
    }
    // optional strength in base class
    if (_e != inert) {
      // need to assign it a vector first!
      std::vector<S> new_s;
      new_s.resize(3*_n);
      for (size_t i=0; i<3*_n; ++i) {
        new_s[i] = (-1.0 + 2.0*(S)rand()/(S)RAND_MAX) / (S)_n;
      }
      this->s = std::move(new_s);
    }
    // velocity in base class
    this->u.resize(3*_n);
    // optional velgrads here
    if (_m == lagrangian) {
      std::vector<S> new_ug;
      new_ug.resize(9*_n);
      ug = std::move(new_ug);
    }
  }

  std::optional<std::vector<S>>& get_velgrad() { return ug; }

  void zero_vels() {
    // must explicitly call the method in the base class
    ElementBase<S>::zero_vels();
    // and specialize
    if (ug) {
      for( S& val : *ug ) val = 0.0;
    }
  }
  void move(const double _dt) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_dt);
    // and specialize
    if (this->M == lagrangian and ug) {
      std::cout << "  Stretching" << to_string() << std::endl;

      // get pointers to the right part of the vectors
      std::vector<S> all_ug = *ug;
      std::vector<S> all_s = *this->s;

      for (size_t i=0; i<this->n; ++i) {
        std::array<S,3> wdu = {0.0};
        S* this_ug = &all_ug[9*i];
        S* this_s = &all_s[3*i];

        // compute stretch term
        // note that multiplying by the transpose may maintain linear impulse better, but
        //   severely underestimates stretch!
        wdu[0] = this_s[0]*this_ug[0] + this_s[1]*this_ug[3] + this_s[2]*this_ug[6];
        wdu[1] = this_s[0]*this_ug[1] + this_s[1]*this_ug[4] + this_s[2]*this_ug[7];
        wdu[2] = this_s[0]*this_ug[2] + this_s[1]*this_ug[5] + this_s[2]*this_ug[8];

        // update elongation
        const S circmag = std::sqrt(this_s[0]*this_s[0] + this_s[1]*this_s[1] + this_s[2]*this_s[2]);
        const S elongfactor = (S)_dt * (this_s[0]*wdu[0] + this_s[1]*wdu[1] + this_s[2]*wdu[2]) / circmag;
        elong[i] *= 1.0 + elongfactor;

        // add Cottet SFS into stretch term (after elongation)

        // update strengths
        this_s[0] += _dt * wdu[0];
        this_s[1] += _dt * wdu[1];
        this_s[2] += _dt * wdu[2];
      }
    }
  }

  std::string to_string() const {
    std::string retstr = ElementBase<S>::to_string() + " Points";
    if (ug) retstr += " with grads";
    return retstr;
  }

protected:
  //Geometry<0,0,S> x;    // num dimensions, num order, storage type
  // movement
  //std::optional<Body&> b;
  // state vector
  std::vector<S> elong;   // scalar elongation
  // time derivative of state vector
  std::optional<std::vector<S>> ug;   // velocity gradients
};


struct SmoothPanelStorage {
};

// 1-D elements
template <class S>
class Panels: public ElementBase<S> {
public:
  Panels(const size_t _n, const elem_t _e, const move_t _m) :
      ElementBase<S>(_n, _e, _m) {
    // this initialization specific to Panels - a circle
    this->x.resize(2*_n);
    this->r.resize(_n);
    for (size_t i=0; i<_n; ++i) {
      this->x[2*i+0] = 0.5 * cos(2.0*i*M_PI/_n);
      this->x[2*i+1] = 0.5 * sin(2.0*i*M_PI/_n);
      this->r[i] = 0.01;
    }
    // initialize indices to nodes
    idx.resize(2*_n);
    for (size_t i=0; i<_n; ++i) {
      idx[2*i] = i;
      idx[2*i+1] = i+1;
    }
    idx[2*_n-1] = 0;
    // just size vels
    this->u.resize(2*_n);
    // and generate panel strengths
    if (_e != inert) {
      // need to assign it a vector first!
      std::vector<S> new_s;
      new_s.resize(_n);
      for (size_t i=0; i<_n; ++i) {
        new_s[i] = 2.0 * sin(2.0*(i+0.5)*M_PI/_n);
      }
      this->s = std::move(new_s);
    }
  }

  const std::vector<uint16_t>& get_idx() const { return idx; }

  std::string to_string() const {
    return ElementBase<S>::to_string() + " Panels";
  }

protected:
  // geometry
  //std::vector<S> x;
  //std::vector<S> r;
  //std::variant<std::vector<uint16_t>, std::vector<uint32_t>> idx;
  std::vector<uint16_t> idx;
  // velocity
  //std::vector<S> u;
  // movement
  //std::optional<Body&> b;
  // curved panels need: normals
  std::optional<std::vector<S>> norm;
};


// velocity influence functions

template <class S, class A>
static inline void kernel_0_0 (const S* __restrict__ sx, const S __restrict__ sr, const S* __restrict__ ss,
                               const S* __restrict__ tx, const S __restrict__ tr, A* __restrict__ tu) {
  // 30 flops
  const S dx = tx[0] - sx[0];
  const S dy = tx[1] - sx[1];
  const S dz = tx[2] - sx[2];
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
  r2 = 1.0 / (r2*std::sqrt(r2));
  const S dxxw = dz*ss[1] - dy*ss[2];
  const S dyxw = dx*ss[2] - dz*ss[0];
  const S dzxw = dy*ss[0] - dx*ss[1];
  tu[0] += r2 * dxxw;
  tu[1] += r2 * dyxw;
  tu[2] += r2 * dxxw;
}

template <class S, class A>
static inline void kernel_0_0g (const S* __restrict__ sx, const S __restrict__ sr, const S* __restrict__ ss,
                                const S* __restrict__ tx, const S __restrict__ tr, A* __restrict__ tu) {
  // 30 flops
  const S dx = tx[0] - sx[0];
  const S dy = tx[1] - sx[1];
  const S dz = tx[2] - sx[2];
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
  r2 = 1.0 / (r2*std::sqrt(r2));
  S dxxw = dz*ss[1] - dy*ss[2];
  S dyxw = dx*ss[2] - dz*ss[0];
  S dzxw = dy*ss[0] - dx*ss[1];
  tu[0] += r2 * dxxw;
  tu[1] += r2 * dyxw;
  tu[2] += r2 * dxxw;

  // HACK - you need to figure out what this term is
  const S bbb = r2 / std::sqrt(r2);
  // continuing with grads - this section is 33 flops
  dxxw *= bbb;
  dyxw *= bbb;
  dzxw *= bbb;
  tu[3] += dx*dxxw;
  tu[4] += dx*dyxw + ss[2]*r2;
  tu[5] += dx*dzxw - ss[1]*r2;
  tu[6] += dy*dxxw - ss[2]*r2;
  tu[7] += dy*dyxw;
  tu[8] += dy*dzxw + ss[0]*r2;
  tu[9] += dz*dxxw + ss[1]*r2;
  tu[10] += dz*dyxw - ss[0]*r2;
  tu[11] += dz*dzxw;
}

template <class S, class A>
void points_affect_points (Points<S> const& src, Points<S>& targ) {
  std::cout << "    0_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::vector<S>& sx = src.get_pos();
  const std::vector<S>& sr = src.get_rad();
  const std::vector<S>& ss = src.get_str();
  const std::vector<S>& tx = targ.get_pos();
  const std::vector<S>& tr = targ.get_rad();
  std::vector<S>& tu = targ.get_vel();
  auto& opttug = targ.get_velgrad();

  // here is where we can dispatch on solver type, grads-or-not, core function, etc.

  // dispatch on presence of val grads
  if (opttug) {
    //std::vector<S>& tug = targ.get_vel();
    // get the pointer from the optional
    auto tug = *opttug;

    // velocity+grads kernel
    for (size_t i=0; i<targ.getn(); ++i) {
      std::array<A,12> accum = {0.0};
      for (size_t j=0; j<src.getn(); ++j) {
        kernel_0_0g<S,A>(&sx[3*j], sr[j], &ss[3*j],
                         &tx[3*i], tr[i], accum.data());
      }
      tu[3*i+0] += accum[0];
      tu[3*i+1] += accum[1];
      tu[3*i+2] += accum[2];
      tug[9*i+0] += accum[3];
      tug[9*i+1] += accum[4];
      tug[9*i+2] += accum[5];
      tug[9*i+3] += accum[6];
      tug[9*i+4] += accum[7];
      tug[9*i+5] += accum[8];
      tug[9*i+6] += accum[9];
      tug[9*i+7] += accum[10];
      tug[9*i+8] += accum[11];
    }

  } else {
    // velocity-only kernel
    for (size_t i=0; i<targ.getn(); ++i) {
      std::array<A,3> accum = {0.0};
      for (size_t j=0; j<src.getn(); ++j) {
        kernel_0_0<S,A>(&sx[3*j], sr[j], &ss[3*j],
                        &tx[3*i], tr[i], accum.data());
      }
      tu[3*i+0] += accum[0];
      tu[3*i+1] += accum[1];
      tu[3*i+2] += accum[2];
    }
  }
}

//
// analytic influence of 2d linear constant-strength vortex panel on target point
//   ignoring the 1/2pi factor, which will be multiplied later
//   40 flops average
//
template <class S, class A>
static inline void kernel_1_0 (const S* __restrict__ sx0, const S* __restrict__ sx1, const S str,
                               const S* __restrict__ tx, A* __restrict__ tu) {

  // side lengths of the triangle s0, s1, t
  const S rij2  = std::pow(tx[0]-sx0[0],2) + std::pow(tx[1]-sx0[1],2);
  const S rij   = std::sqrt(rij2);
  const S rij12 = std::pow(tx[0]-sx1[0],2) + std::pow(tx[1]-sx1[1],2);
  const S rij1  = std::sqrt(rij12);
  //std::cout << "rij is " << rij << " and rijp1 is " << rij1 << std::endl;
  const S vstar = std::log(rij/rij1);
  S ustar = std::atan2(tx[0]-sx1[0], tx[1]-sx1[1]) - std::atan2(tx[0]-sx0[0], tx[1]-sx0[1]);
  //std::cout << "ustar started off as " << ustar << std::endl;
  if (ustar < -M_PI) ustar += 2.*M_PI;
  if (ustar > M_PI) ustar -= 2.*M_PI;
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  const S px    = sx1[0]-sx0[0];
  const S py    = sx1[1]-sx0[1];
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  const S mult  = str / std::sqrt(std::pow(px,2) + std::pow(py,2));
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and multiply by vortex sheet strength
  tu[0] += mult*velx;
  tu[1] += mult*vely;
}

template <class S, class A>
void panels_affect_points (Panels<S> const& src, Points<S>& targ) {
  std::cout << "    1_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::vector<S>& sx = src.get_pos();
  //const std::vector<S>& sr = src.get_rad();
  const std::vector<uint16_t>& si = src.get_idx();
  const std::vector<S>& ss = src.get_str();
  const std::vector<S>& tx = targ.get_pos();
  //const std::vector<S>& tr = targ.get_rad();
  std::vector<S>& tu = targ.get_vel();

  for (size_t i=0; i<targ.getn(); ++i) {
    std::array<A,2> accum = {0.0};
    for (size_t j=0; j<src.getn(); ++j) {
      kernel_1_0<S,A>(&sx[2*si[2*j]], &sx[2*si[2*j+1]], ss[j],
                      &tx[2*i], accum.data());
    }
    tu[2*i]   += accum[0];
    tu[2*i+1] += accum[1];
  }
}

template <class S, class A>
void points_affect_panels (Points<S> const& src, Panels<S>& targ) {
  std::cout << "    0_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::vector<S>& sx = src.get_pos();
  const std::vector<S>& ss = src.get_str();
  const std::vector<S>& tx = targ.get_pos();
  const std::vector<uint16_t>& ti = targ.get_idx();
  std::vector<S>& tu = targ.get_vel();

  for (size_t i=0; i<targ.getn(); ++i) {
    std::array<A,2> accum = {0.0};
    for (size_t j=0; j<src.getn(); ++j) {
      // note that this is the same kernel as panels_affect_points!
      kernel_1_0<S,A>(&tx[2*ti[2*i]], &tx[2*ti[2*i+1]], ss[j],
                      &sx[2*j], accum.data());
    }
    // we use it backwards, so the resulting velocities are negative
    tu[2*i]   -= accum[0];
    tu[2*i+1] -= accum[1];
  }
}

template <class S, class A>
void panels_affect_panels (Panels<S> const& src, Panels<S>& targ) {
  std::cout << "    1_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  // not sure how to do this - find field points of one and apply a function above?
}


// helper struct for dispatching through a variant
template <class A>
struct InfluenceVisitor {
  // source collection, target collection
  void operator()(Points<float> const& src, Points<float>& targ) { points_affect_points<float,A>(src, targ); } 
  void operator()(Panels<float> const& src, Points<float>& targ) { panels_affect_points<float,A>(src, targ); } 
  void operator()(Points<float> const& src, Panels<float>& targ) { points_affect_panels<float,A>(src, targ); } 
  void operator()(Panels<float> const& src, Panels<float>& targ) { panels_affect_panels<float,A>(src, targ); } 
};


// helper alias for any type of collection of elements
using Collection = std::variant<Points<float>,
				Panels<float>>;


// helper function
std::string to_string(Collection& _c) {
  return std::visit([=](auto& elem) { return elem.to_string(); }, _c);
}





// execution starts here

int main(int argc, char const *argv[]) {

  const solution_t solver = direct_cpu;
  std::array<double,Dimensions> fs = {0.0, 0.0, 0.0};
  double time = 0.0;
  const double dt = 0.01;
 
  std::vector<Collection> vort;		// the free vorticity
  std::vector<Collection> bdry;		// all boundaries
  std::vector<Collection> fldpt;	// tracers and field points
 
  vort.push_back(Points<float>(10000, active, lagrangian));	// vortons
  //fldpt.push_back(Points<float>(2000, inert, lagrangian));	// tracer particles
  //fldpt.push_back(Points<float>(100, inert, fixed));		// static field points
  //bdry.push_back(Panels<float>(500, reactive, bodybound));	// panels

  // need this for dispatching velocity influence calls, template param is accumulator type
  // should the solution_t be an argument to the constructor?
  InfluenceVisitor<double> visitor;


  // one-half diffusion step
  //tbd


  // this is one Euler convection step - how would we do one 2nd order step?
  // can we have an object for different forward integrators?

  if (bdry.size() > 0) std::cout << std::endl << "Solving for BEM RHS" << std::endl << std::endl;

  // loop over boundary collections
  for (auto &targ: bdry) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl;
    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
    // accumulate from vorticity
    for (auto &src: vort) {
      std::visit(visitor, src, targ);
    }
  }
 
  if (bdry.size() > 0) std::cout << std::endl << "Solving for BEM matrix" << std::endl << std::endl;

  // loop over boundary collections
  for (auto &targ: bdry) {
    std::cout << "  Solving for influence coefficients on" << to_string(targ) << std::endl;
    // assemble from all boundaries
    for (auto &src: bdry) {
      std::visit(visitor, src, targ);
    }
  }
 
  //std::cout << std::endl << "Solving BEM for strengths" << std::endl << std::endl;

  if (vort.size()+fldpt.size() > 0) std::cout << std::endl << "Solving for velocities" << std::endl << std::endl;

  // find the influence on every vorticity element
  for (auto &targ: vort) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl << std::flush;
    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
    // accumulate from vorticity
    for (auto &src: vort) {
      // how do I specify the solver?
      std::visit(visitor, src, targ);
    }
    // accumulate from boundaries
    for (auto &src: bdry) {
      std::visit(visitor, src, targ);
    }
    // add freestream and divide by 2pi
    std::visit([=](auto& elem) { elem.finalize_vels(fs); }, targ);
  }
 
  // find the influence on every field point/tracer element
  for (auto &targ: fldpt) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl;
    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
    // accumulate from vorticity
    for (auto &src: vort) {
      std::visit(visitor, src, targ);
    }
    // accumulate from boundaries
    for (auto &src: bdry) {
      std::visit(visitor, src, targ);
    }
    // add freestream and divide by 2pi
    std::visit([=](auto& elem) { elem.finalize_vels(fs); }, targ);
  }

  std::cout << std::endl << "Convection step" << std::endl << std::endl;

  // move every movable element
  for (auto &coll: vort) {
    std::visit([=](auto& elem) { elem.move(dt); }, coll);
  }
  for (auto &coll: bdry) {
    std::visit([=](auto& elem) { elem.move(dt); }, coll);
  }
  for (auto &coll: fldpt) {
    std::visit([=](auto& elem) { elem.move(dt); }, coll);
  }
 
  // one-half diffusion step
  //tbd

  // write output, restart files
  time += dt;

  std::cout << std::endl << "Done" << std::endl << std::endl;

  return 0;
}

