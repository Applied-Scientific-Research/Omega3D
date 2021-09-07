/*
 * Simulation.cpp - a class to control a 3D vortex particle sim
 *
 * (c)2017-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#include "Simulation.h"
#include "Reflect.h"
#include "BEMHelper.h"
#include "VtkXmlHelper.h"
#include "Split.h"
#include "GuiHelper.h"

#include <cassert>
#include <cmath>
#include <cfenv> // Catch fp exceptions
#include <limits>
#include <variant>

#ifdef _WIN32
#pragma STDC FENV_ACCESS ON // For fp exceptions
#endif


// constructor
Simulation::Simulation()
  : re(100.0),
    dt(0.01),
    fs{0.0,0.0,0.0},
    bodies(),
    vort(),
    bdry(),
    fldpt(),
    bem(),
    diff(),
    conv(),
    sf(),
    description(),
    time(0.0),
    output_dt(0.0),
    end_time(100.0),
    use_end_time(false),
    overlap_ratio(1.5),
    core_size_ratio(std::sqrt(8.0)),
    nstep(0),
    use_max_steps(false),
    max_steps(100),
    auto_start(false),
    quit_on_stop(false),
    sim_is_initialized(false),
    step_has_started(false),
    step_is_finished(false)
  {}

// addresses for use in imgui
float* Simulation::addr_re() { return &re; }
float* Simulation::addr_dt() { return &dt; }
float* Simulation::addr_fs() { return fs; }

// getters
float Simulation::get_re() const { return re; }
float Simulation::get_dt() const { return dt; }
float Simulation::get_hnu() const { return std::sqrt(dt/re); }
float Simulation::get_ips() const { return core_size_ratio * get_hnu(); }
float Simulation::get_vdelta() const { return overlap_ratio * get_ips(); }
float Simulation::get_time() const { return (float)time; }
float Simulation::get_end_time() const { return (float)end_time; }
bool Simulation::using_end_time() const { return use_end_time; }
size_t Simulation::get_nstep() const { return nstep; }
size_t Simulation::get_max_steps() const { return max_steps; }
bool Simulation::using_max_steps() const { return use_max_steps; }
float Simulation::get_output_dt() const { return (float)output_dt; }
std::string Simulation::get_description() { return description; }
bool Simulation::autostart() { return auto_start; }
bool Simulation::quitonstop() { return quit_on_stop; }

// setters
void Simulation::set_description(const std::string desc) { description = desc; }
void Simulation::set_end_time(const double net) { end_time = net; use_end_time = true; }
void Simulation::unset_end_time() { use_end_time = false; }
void Simulation::set_max_steps(const size_t nms) { max_steps = nms; use_max_steps = true; }
void Simulation::unset_max_steps() { use_max_steps = false; }
void Simulation::set_output_dt(const double nodt) { output_dt = nodt; }
void Simulation::set_auto_start(const bool autos) { auto_start = autos; }
void Simulation::set_quit_on_stop(const bool qos) { quit_on_stop = qos; }

// access status file
void Simulation::set_status_file_name(const std::string _fn) { sf.set_filename(_fn); }
std::string Simulation::get_status_file_name() { return sf.get_filename(); }

// status
size_t Simulation::get_npanels() {
  size_t n = 0;
  for (auto &coll: bdry) {
    //std::visit([&n](auto& elem) { n += elem.get_npanels(); }, coll);
    // only proceed if the last collection is Surfaces
    if (std::holds_alternative<Surfaces<float>>(coll)) {
      Surfaces<float>& surf = std::get<Surfaces<float>>(coll);
      n += surf.get_npanels();
    }
  }
  return n;
}

size_t Simulation::get_nparts() {
  size_t n = 0;
  for (auto &coll: vort) {
    std::visit([&n](auto& elem) { n += elem.get_n(); }, coll);
  }
  return n;
}

size_t Simulation::get_nfldpts() {
  size_t n = 0;
  for (auto &coll : fldpt) {
    std::visit([&n](auto& elem) { n += elem.get_n(); }, coll);
  }
  return n;
}

// like a setter
void Simulation::set_re_for_ips(const float _ips) {
  re = std::pow(core_size_ratio, 2) * dt / pow(_ips, 2);
  diff.set_diffuse(false);
}

void Simulation::set_diffuse(const bool _do_diffuse) {
  diff.set_diffuse(_do_diffuse);
}

void Simulation::set_amr(const bool _do_amr) {
  diff.set_amr(_do_amr);
}


//
// json read/write
//

// read "simparams" json object
void
Simulation::from_json(const nlohmann::json j) {

  if (j.find("nominalDt") != j.end()) {
    dt = j["nominalDt"];
    std::cout << "  setting dt= " << dt << std::endl;
  }

  if (j.find("outputDt") != j.end()) {
    output_dt = j["outputDt"];
    std::cout << "  setting output dt= " << output_dt << std::endl;
  }

  //if (j.find("nominalDx") != j.end()) {
  //  dx = j["nominalDx"];
  //  std::cout << "  setting dx= " << dx << std::endl;
  //}

  if (j.find("maxSteps") != j.end()) {
    use_max_steps = true;
    max_steps = j["maxSteps"];
    std::cout << "  setting max_steps= " << max_steps << std::endl;
  } else {
    use_max_steps = false;
  }

  if (j.find("endTime") != j.end()) {
    use_end_time = true;
    end_time = j["endTime"];
    std::cout << "  setting end_time= " << end_time << std::endl;
  } else {
    use_end_time = false;
  }

  if (j.find("overlapRatio") != j.end()) {
    overlap_ratio = j["overlapRatio"];
    std::cout << "  setting overlap ratio= " << overlap_ratio << std::endl;
  }

  if (j.find("coreSizeRatioSqrd") != j.end()) {
    core_size_ratio = std::sqrt((float)j["coreSizeRatioSqrd"]);
    std::cout << "  setting core size ratio (nominal separation over h_nu) = " << core_size_ratio << std::endl;
  }

  // Convection will find and set "timeOrder"
  conv.from_json(j);

  // Diffusion will find and set "viscous", "VRM" and "AMR" parameters
  diff.from_json(j);
}

// create and write a json object for "simparams"
nlohmann::json
Simulation::to_json() const {
  nlohmann::json j;

  j["nominalDt"] = dt;
  j["outputDt"] = output_dt;
  if (using_max_steps()) j["maxSteps"] = get_max_steps();
  if (using_end_time()) j["endTime"] = get_end_time();
  j["overlapRatio"] = overlap_ratio;
  j["coreSizeRatioSqrd"] = std::pow(core_size_ratio,2);

  // Convection will write "timeOrder"
  conv.add_to_json(j);

  // Diffusion will write "viscous", "VRM" and "AMR" parameters
  diff.add_to_json(j);

  return j;
}

// set "flowparams" json object
void
Simulation::flow_from_json(const nlohmann::json j) {

  if (j.find("Re") != j.end()) {
    re = j["Re"];
    std::cout << "  setting re= " << re << std::endl;
  }
  if (j.find("Uinf") != j.end()) {
    // eventually support an expression for Uinf instead of just a single float
    std::vector<float> new_fs = {0.0, 0.0, 0.0};
    new_fs.resize(Dimensions);
    if (j["Uinf"].is_array()) {
      new_fs = j["Uinf"].get<std::vector<float>>();
    } else if (j["Uinf"].is_number()) {
      new_fs[0] = j["Uinf"].get<float>();
    }
    for (size_t i=0; i<Dimensions; ++i) fs[i] = new_fs[i];
    std::cout << "  setting freestream to " << fs[0] << " " << fs[1] << " " << fs[2] << std::endl;
  }
}

// create and write a json object for "flowparams"
nlohmann::json
Simulation::flow_to_json() const {
  nlohmann::json j;

  j["Re"] = re;
  j["Uinf"] = {fs[0], fs[1], fs[2]};

  return j;
}

// set "runtime" json object
void
Simulation::runtime_from_json(const nlohmann::json j) {

  sf.from_json(j);

  if (j.find("autoStart") != j.end()) {
    bool autostart = j["autoStart"];
    set_auto_start(autostart);
    std::cout << "  autostart? " << autostart << std::endl;
  }
  if (j.find("quitOnStop") != j.end()) {
    bool qos = j["quitOnStop"];
    set_quit_on_stop(qos);
    std::cout << "  quit on stop? " << qos << std::endl;
  }
}

// create and write a json object for "runtime"
nlohmann::json
Simulation::runtime_to_json() const {
  nlohmann::json j;

  sf.add_to_json(j);

  j["autoStart"] = auto_start;
  j["quitOnStop"] = quit_on_stop;

  return j;
}


#ifdef USE_IMGUI
//
// ImGui code
//
void Simulation::draw_advanced() {

  // set the execution environment in Convection.h
  conv.draw_advanced();

  // set the diffusion parameters in Diffusion.h
  diff.draw_advanced();
}
#endif


#ifdef USE_GL
//
// OpenGL-specific code
//

void Simulation::updateGL() {
  for (auto &coll : vort) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
  for (auto &coll : bdry) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
  for (auto &coll : fldpt) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
}

void Simulation::drawGL(std::vector<float>& _mviewmat,
                        std::vector<float>& _projmat,
                        RenderParams&       _rparams) {

  if (step_is_finished) {
    _rparams.tracer_size = get_ips() * _rparams.tracer_scale;
    for (auto &coll : vort) {
      std::visit([&](auto& elem) { elem.drawGL(_mviewmat, _projmat, _rparams, get_vdelta()); }, coll);
    }
    for (auto &coll : bdry) {
      std::visit([&](auto& elem) { elem.drawGL(_mviewmat, _projmat, _rparams, get_vdelta()); }, coll);
    }
    for (auto &coll : fldpt) {
      std::visit([&](auto& elem) { elem.drawGL(_mviewmat, _projmat, _rparams, get_vdelta()); }, coll);
    }
  }
}
#endif

#ifdef USE_OGL_COMPUTE
//
// Compute shader code
//

void Simulation::initGLcs() {

  // generate the opengl state object with space for 6 vbos and 2 shader programs
  cgl = std::make_shared<GlComputeState<STORE>>(6,2);

  // Load and create the blob-drawing shader program
  cgl->spo[0] = create_ptptvelgrad_program();

  // do the same for the zeroing program
  cgl->spo[1] = create_initvelgrad_program();

  // Identify and bind
  for (GLuint i=0; i<6; ++i) {
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[i]);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, i, cgl->vbo[i]);
  }

  // locate where the offsets and counts go
  cgl->source_offset_attr = glGetUniformLocation(cgl->spo[0], "isrc");
  cgl->source_count_attr  = glGetUniformLocation(cgl->spo[0], "nsrc");
  cgl->target_offset_attr = glGetUniformLocation(cgl->spo[0], "itarg");
  cgl->target_count_attr  = glGetUniformLocation(cgl->spo[0], "ntarg");

  cgl->target_offset_attr_z = glGetUniformLocation(cgl->spo[1], "itarg");
  cgl->target_count_attr_z  = glGetUniformLocation(cgl->spo[1], "ntarg");
}

// upload the new data arrays to the GPU
void Simulation::updateGLcs() {

  assert(cgl && "No GlComputeState object has been made");
  assert(glIsVertexArray(cgl->vao) != GL_FALSE && "GlComputeState VAO is not a VAO");

  cgl->nsrc  = cgl->hsx.size()/4;
  cgl->ntarg = cgl->htx.size()/4;
  const GLsizeiptr szsrc = cgl->hsx.size() * sizeof(STORE);
  const GLsizeiptr sztrg = cgl->htx.size() * sizeof(STORE);

  if (cgl->ntarg > 0) {
    glBindVertexArray(cgl->vao);

    // set the zero-vels shader and run it
    //glUseProgram(cgl->spo[1]);
    //glUniform1ui(cgl->target_offset_attr_z, (const GLuint)0);
    //glUniform1ui(cgl->target_count_attr_z,  (const GLuint)cgl->ntarg);

    // allocate space for the results
    //glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[3]);
    //glBufferData(GL_SHADER_STORAGE_BUFFER, sztrg, nullptr, GL_DYNAMIC_DRAW);
    //glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[4]);
    //glBufferData(GL_SHADER_STORAGE_BUFFER, sztrg, nullptr, GL_DYNAMIC_DRAW);
    //glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[5]);
    //glBufferData(GL_SHADER_STORAGE_BUFFER, sztrg, nullptr, GL_DYNAMIC_DRAW);

    // do the work
    //GLuint group_count = ((GLuint)cgl->ntarg + 128 - 1) / 128;
    //glDispatchCompute( group_count, 1, 1 );
    //glMemoryBarrier( GL_SHADER_STORAGE_BARRIER_BIT );


    // and do the same for the vel evaluation shader
    glUseProgram(cgl->spo[0]);

    // send the current values
    glUniform1ui(cgl->source_offset_attr, (const GLuint)0);
    glUniform1ui(cgl->source_count_attr,  (const GLuint)cgl->nsrc);
    glUniform1ui(cgl->target_offset_attr, (const GLuint)0);
    glUniform1ui(cgl->target_count_attr,  (const GLuint)cgl->ntarg);

    // send the arrays
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[0]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, szsrc, cgl->hsx.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[1]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, szsrc, cgl->hss.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[2]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sztrg, cgl->htx.data(), GL_DYNAMIC_DRAW);

    // and allocate space for the results
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[3]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sztrg, nullptr, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[4]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sztrg, nullptr, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[5]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sztrg, nullptr, GL_DYNAMIC_DRAW);

    glBindVertexArray(0);
  }
}

// download the result data arrays from the GPU
void Simulation::retrieveGLcs() {

  // bind and prepare
  glBindVertexArray(cgl->vao);
  glUseProgram(cgl->spo[0]);

  // get the three results arrays
  const GLsizeiptr szres = cgl->hr1.size() * sizeof(STORE);
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[3]);
  glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, szres, cgl->hr1.data());
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[4]);
  glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, szres, cgl->hr2.data());
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, cgl->vbo[5]);
  glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, szres, cgl->hr3.data());

  glBindVertexArray(0);
}

void Simulation::computeGL() {

  assert(cgl && "No GlComputeState object has been made");

  // check for access to lock, if access, continue
  if (cgl->cstate.load() == no_compute) return;

  // if computation is ready and not running, get it ready
  if (cgl->cstate.load() == begin_compute) {
    // upload the compacted arrays to the GPU
    updateGLcs();
    cgl->cstate.store(computing);
  }

  if (cgl->cstate.load() == computing) {
    // for now, sliver (block) only over targets
    static GLuint sliver_index = 0;

    if (sliver_index == 0) std::cout << "    running sliver." << std::flush;
    else std::cout << "." << std::flush;

    // estimate total work and number of slivers (2 TFlop/s and 60 fps)
    float total_work = 1.0 + (float)cgl->ntarg * (float)cgl->nsrc * (float)flops_0v_0bg<float>();
    //GLuint num_slivers = 1 + (GLuint)(total_work / (150.e+9/30.0));
    GLuint num_slivers = 1 + (GLuint)(total_work / (2.e+12/30.0));
    //GLuint num_slivers = 3;

    // when slivering over targets...
    // calculate work group count
    GLuint total_targ_grps = ((GLuint)cgl->ntarg + 128 - 1) / 128;
    //std::cout << "  total_targ_grps is " << total_targ_grps << std::endl;
    //GLuint groups_per_sliver = (total_targ_grps + num_slivers - 1) / num_slivers;
    // find offset and count for this sliver
    //GLuint group_offset = sliver_index*groups_per_sliver;
    //GLuint targ_offset = 128*group_offset;
    //GLuint group_count = std::min(groups_per_sliver, total_targ_grps-group_offset);
    //GLuint targ_count = std::min(128*group_count, cgl->ntarg-targ_offset);
    //std::cout << "    sliver of " << group_count << " target groups offset by " << group_offset << std::endl << std::flush;
    //std::cout << "           or " << cgl->nsrc << " sources on " << targ_count << " targets" << std::endl << std::flush;

    // when slivering over sources...
    GLuint total_src_grps = ((GLuint)cgl->nsrc + 128 - 1) / 128;
    GLuint groups_per_sliver = (total_src_grps + num_slivers - 1) / num_slivers;
    // find offset and count for this sliver
    GLuint group_offset = sliver_index*groups_per_sliver;
    GLuint src_offset = 128*group_offset;
    GLuint group_count = std::min(groups_per_sliver, total_src_grps-group_offset);
    GLuint src_count = std::min(128*group_count, cgl->nsrc-src_offset);

    // set up the compute shader program
    glBindVertexArray(cgl->vao);
    glUseProgram(cgl->spo[0]);

    // sliver over targets, each target integrates over all sources
    //glUniform1ui(cgl->source_offset_attr, (const GLuint)0);
    //glUniform1ui(cgl->source_count_attr,  (const GLuint)cgl->nsrc);
    //glUniform1ui(cgl->target_offset_attr, (const GLuint)targ_offset);
    //glUniform1ui(cgl->target_count_attr,  (const GLuint)targ_count);

    // sliver over sources, every target adds some influence
    glUniform1ui(cgl->source_offset_attr, (const GLuint)src_offset);
    glUniform1ui(cgl->source_count_attr,  (const GLuint)src_count);
    glUniform1ui(cgl->target_offset_attr, (const GLuint)0);
    glUniform1ui(cgl->target_count_attr,  (const GLuint)cgl->ntarg);
    //std::cout << "    sliver of sources " << src_offset << " " << src_count << " on targets " << 0 << " " << cgl->ntarg << std::endl << std::flush;

    // do the work
    glDispatchCompute( total_targ_grps, 1, 1 );
    glMemoryBarrier( GL_SHADER_STORAGE_BARRIER_BIT );
    //glFinish();	// do we need this? probably not
    glBindVertexArray(0);

    // increment the offset
    sliver_index++;

    // only trigger the atomic once all slivers are complete
    if (sliver_index == num_slivers) {
      // reset the offset for the next time
      sliver_index = 0;
      std::cout << std::endl << std::flush;

      // atomically indicate completion
      cgl->cstate.store(compute_done);
    }
  }

  if (cgl->cstate.load() == compute_done) {
    // download the compacted arrays from the GPU
    retrieveGLcs();
    //glFinish();	// do we need this? probably not

    // tell flow objects to update their values to the GPU
    updateGL();

    // tell the other thread that we're done
    cgl->cstate.store(no_compute);
  }
}
#endif

bool Simulation::is_initialized() { return sim_is_initialized; }

void Simulation::set_initialized() { sim_is_initialized = true; }

void Simulation::reset() {

  // must wait for step() to complete, if it's still working
  if (stepfuture.valid()) {
    stepfuture.wait();
    stepfuture.get();
  }

  // now reset everything else
  time = 0.0;
  nstep = 0;
  vort.clear();
  bdry.clear();
  fldpt.clear();
  bem.reset();
  sf.reset_sim();
  sim_is_initialized = false;
  step_has_started = false;
  step_is_finished = false;
}

void Simulation::clear_bodies() {
  bodies.clear();
}

// Write a set of vtu files for the particles and panels
std::vector<std::string> Simulation::write_vtk(const int _index,
                                               const bool _do_bdry,
                                               const bool _do_flow,
                                               const bool _do_measure) {

  std::cout << "Inside Simulation::write_vtk at t=" << time << std::endl;

  // solve the BEM (before any VTK or status file output)
  //std::cout << "Updating element vels" << std::endl;
#ifndef USE_OGL_COMPUTE
  std::array<double,3> thisfs = {fs[0], fs[1], fs[2]};
  //clear_inner_layer<STORE>(1, bdry, vort, 1.0/std::sqrt(2.0*M_PI), get_ips());
  solve_bem<STORE,ACCUM,Int>(time, thisfs, vort, bdry, bem);

  if (_do_flow)    conv.find_vels(thisfs, vort, bdry, vort, true);
  if (_do_measure) conv.find_vels(thisfs, vort, bdry, fldpt, true);
  if (_do_bdry)    conv.find_vels(thisfs, vort, bdry, bdry, true);
#endif

  // may eventually want to avoid clobbering by maintaining an internal count of the
  //   number of simulations run from this execution of the GUI
  std::vector<std::string> files;

  // if a positive number was passed in, use that instead of the current step
  size_t stepnum = 0;
  if (_index < 0) {
    stepnum = nstep;
  } else {
    stepnum = (size_t)_index;
  }

  // ask Vtk to write files for each collection
  if (_do_flow)    write_vtk_files<float>(vort, stepnum, time, files);
  if (_do_measure) write_vtk_files<float>(fldpt, stepnum, time, files);
  if (_do_bdry)    write_vtk_files<float>(bdry, stepnum, time, files);

  return files;
}

//
// Check all aspects of the initialization for conditions that prevent a run from starting
//
std::string Simulation::check_initialization() {
  std::string retstr;

  // are any flow features particle generators? - HOW DO WE DO THIS?
  const bool are_generators = false;

  // Check for no bodies and no particles
  if (get_npanels() == 0 and get_nparts() == 0 and not are_generators) {
    retstr.append("No flow features and no bodies. Add one or both, reset, and run.\n");
  }

  // Check for a body and no particles
  if (get_npanels() > 0 and get_nparts() == 0) {

    const bool zero_freestream = (fs[0]*fs[0]+fs[1]*fs[1] < std::numeric_limits<float>::epsilon());
    const bool no_body_movement = not do_any_bodies_move();
    const bool all_zero_bcs = not any_nonzero_bcs();

    // AND no freestream
    if (zero_freestream and no_body_movement and all_zero_bcs) {
      retstr.append("No flow features, zero freestream speed, no movement, and no driven boundaries - try adding one of these.\n");
      return retstr;
    }

    // AND no viscosity
    if (not diff.get_diffuse()) {
      retstr.append("You have a solid body, but no diffusion. It will not shed vorticity. Turn on viscosity or add a flow feature, reset, and run.\n");
    }
  }

  // Check for very large BEM problem
  if (get_npanels() > 21000) {
    retstr.append("Boundary features have too many panels, program will run out of memory. Reduce Reynolds number or increase time step or both.\n");
  }

  // Check for excessive elongation
  float max_elong = 0.0;
  for (auto &coll: vort) {
    // only check particles ("Points")
    if (std::holds_alternative<Points<float>>(coll)) {
      Points<float>& pts = std::get<Points<float>>(coll);
      max_elong = std::max(max_elong, pts.get_max_elong());
    }
    //std::visit([&](auto& elem) { max_elong = std::max(max_elong, elem.get_max_elong()); }, coll);
  }
  if (max_elong > 1.5) retstr.append("Elongation threshold exceeded! Reset and reduce the time step size.\n");

  // Check for vorticity-based Courant number - too high and we're probably messing up
  // later

  return retstr;
}

//
// Check dynamic aspects of the simulation for conditions that should stop the run
//
std::string Simulation::check_simulation() {
  std::string retstr;

  // Are there any dynamic problems in 2D that could blow a run?

  return retstr;
}

//
// check all bodies for movement
//
bool Simulation::do_any_bodies_move() {
  bool some_move = false;
  for (size_t i=0; i<bodies.size(); ++i) {
    auto thisvel = bodies[i]->get_vel(time);
    auto nextvel = bodies[i]->get_vel(time+dt);
    auto thisrot = bodies[i]->get_rotvel_vec(time);
    auto nextrot = bodies[i]->get_rotvel_vec(time+dt);
    if (std::abs(thisvel[0]) + std::abs(thisvel[1]) + std::abs(thisvel[2]) +
        std::abs(thisrot[0]) + std::abs(thisrot[1]) + std::abs(thisrot[2]) +
        std::abs(nextvel[0]) + std::abs(nextvel[1]) + std::abs(nextvel[2]) +
        std::abs(nextrot[0]) + std::abs(nextrot[1]) + std::abs(nextrot[2]) >
        std::numeric_limits<float>::epsilon()) {
      some_move = true;
    }
  }
  return some_move;
}

//
// check all bodies for non-zero BCs
//
bool Simulation::any_nonzero_bcs() {
  bool all_are_zero = true;
  // loop over bdry Collections
  for (auto &src : bdry) {
    float max_bc = std::visit([=](auto& elem) { return elem.get_max_bc_value(); }, src);
    if (std::abs(max_bc) > std::numeric_limits<float>::epsilon()) all_are_zero = false;
  }
  return not all_are_zero;
}

//
// query and get() the future if possible
//
bool Simulation::test_for_new_results() {

  if (not step_has_started) {
    // if we haven't made an async call yet
    return true;

  } else if (is_future_ready(stepfuture)) {
    // if we did, and it's ready to give us the results
    stepfuture.get();

#ifdef USE_GL
    // tell flow objects to update their values to the GPU
    updateGL();
#endif

    // set flag indicating that at least one step has been solved
    step_is_finished = true;
    step_has_started = false;

    return true;
  }

  // async call is not finished, do not try calling it again
  return false;
}

//
// call this from a real-time GUI - will only start the first step if no other steps are being worked on
//
void Simulation::async_first_step() {
  step_has_started = true;
  stepfuture = std::async(std::launch::async, [this](){first_step();});
}

//
// initialize the system so we can start drawing things
//
void Simulation::first_step() {
  std::cout << std::endl << "Taking step " << nstep << " at t=" << time << std::endl;

  // we wind up using this a lot
  std::array<double,3> thisfs = {fs[0], fs[1], fs[2]};

  // this is the first step, just solve BEM and return - it's time=0

  // update BEM and find vels on any particles but DO NOT ADVECT
  conv.advect_1st(time, 0.0, thisfs, get_ips(), vort, bdry, fldpt, bem);

  // and write status file
  dump_stats_to_status();
}

//
// call this from a real-time GUI - will only run a new step if the last one is done
//
void Simulation::async_step() {
  step_has_started = true;
  stepfuture = std::async(std::launch::async, [this](){step();});
}

//
// here's the vortex method: convection and diffusion with operator splitting
//
void Simulation::step() {
  std::cout << std::endl << "Taking step " << nstep << " at t=" << time << " with n=" << get_nparts() << std::endl;

  const bool use_2nd_order_operator_splitting = true;

  // we wind up using this a lot
  std::array<double,3> thisfs = {fs[0], fs[1], fs[2]};

  if (use_2nd_order_operator_splitting) {
    // operator splitting requires one half-step diffuse (use coefficients from previous step, if available)
    diff.step(time, 0.5*dt, re, overlap_ratio, get_vdelta(), thisfs, vort, bdry, bem);
  } else {
    // for simplicity's sake, just run one full diffusion step here
    diff.step(time, dt, re, overlap_ratio, get_vdelta(), thisfs, vort, bdry, bem);
  }

  // advect with no diffusion (must update BEM strengths)
  //conv.advect_1st(time, dt, thisfs, get_ips(), vort, bdry, fldpt, bem);
  conv.advect_2nd(time, dt, thisfs, get_ips(), vort, bdry, fldpt, bem);

  if (use_2nd_order_operator_splitting) {
    // operator splitting requires another half-step diffuse (must compute new coefficients)
    diff.step(time+dt, 0.5*dt, re, overlap_ratio, get_vdelta(), thisfs, vort, bdry, bem);
  }

  // step complete, now split any elongated particles (move this to Split)
  for (auto &coll: vort) {
  
    // but only check particles ("Points")
    if (std::holds_alternative<Points<float>>(coll)) {

      Points<float>& pts = std::get<Points<float>>(coll);
      //std::cout << "    check split for " << pts.get_n() << " particles" << std::endl;
      //std::cout << std::endl;

      // none of these are passed as const, because both may be extended with new particles
      std::array<Vector<float>,Dimensions>&       x = pts.get_pos();
      Vector<float>&                              r = pts.get_rad();
      Vector<float>&                              elong = pts.get_elong();
      std::array<Vector<float>, numStrenPerNode>& s = pts.get_str();

      // last two arguments are: relative distance, allow variable core radii
      (void)split_elongated<float>(x[0], x[1], x[2], r, elong, s[0], s[1], s[2],
                                   diff.get_core_func(),
                                   overlap_ratio,
                                   1.2);

      // we probably have a different number of particles now, resize the u, ug, elong arrays
      pts.resize(r.size());
    }
  }

  // update time
  time += (double)dt;

  // push field points out of objects every few steps
  if (nstep%5 == 0) clear_inner_layer<STORE>(1, bdry, fldpt, (STORE)0.0, (STORE)(0.5*get_ips()));

  // only increment step here!
  nstep++;

  // and write status file
  dump_stats_to_status();
}

//
// close out the step with some work and output to the status file
//
void Simulation::dump_stats_to_status() {
  if (sf.is_active()) {
    // the basics
    sf.append_value("time",(float)time);
    sf.append_value("Nv",(int)get_nparts());

    // more advanced info

    std::array<double,3> thisfs = {fs[0], fs[1], fs[2]};

    // push away particles inside or too close to the body
    //clear_inner_layer<STORE>(1, bdry, vort, 1.0/std::sqrt(2.0*M_PI), get_ips());
    // solve the BEM (before any VTK or status file output)
    solve_bem<STORE,ACCUM,Int>(time, thisfs, vort, bdry, bem);

    // but do we really need to do these?
    //conv.find_vels(thisfs, vort, bdry, vort);
    //conv.find_vels(thisfs, vort, bdry, fldpt);
    //conv.find_vels(thisfs, vort, bdry, bdry);

    // add up the total circulation
    std::array<float,3> tot_circ = {0.0};
    for (auto &src : vort) {
      auto this_circ = std::visit([=](auto& elem) { return elem.get_total_circ(time); }, src);
      for (size_t i=0; i<3; ++i) tot_circ[i] += this_circ[i];
    }
    // then add up the circulation in bodies - DO WE NEED TO RE-SOLVE BEM FIRST?
    for (auto &src : bdry) {
      auto this_circ = std::visit([=](auto& elem) { return elem.get_total_circ(time); }, src);
      for (size_t i=0; i<3; ++i) tot_circ[i] += this_circ[i];
      this_circ = std::visit([=](auto& elem) { return elem.get_body_circ(time); }, src);
      for (size_t i=0; i<3; ++i) tot_circ[i] += this_circ[i];
    }
    sf.append_value("gx",tot_circ[0]);
    sf.append_value("gy",tot_circ[1]);
    sf.append_value("gz",tot_circ[2]);

    // now forces
    std::array<float,Dimensions> impulse = calculate_simple_forces();
    sf.append_value("fx",impulse[0]);
    sf.append_value("fy",impulse[1]);
    if (Dimensions > 2) sf.append_value("fz",impulse[3]);

    // write here
    sf.write_line();
  }
}

// Use impulse method to calculate total forces
std::array<float,Dimensions>
Simulation::calculate_simple_forces() {

  static double last_time = 0.0;
  static std::array<float,Dimensions> last_impulse = {0.0};

  // reset the "last" values if time is zero
  if (time < 0.1*dt) {
    last_time = -dt;
    last_impulse.fill(0.0);
  }

  std::array<float,Dimensions> this_impulse = calculate_total_impulse();

  // find the time derivative of the impulses
  std::array<float,Dimensions> forces;
  for (size_t i=0; i<Dimensions; ++i) forces[i] = (this_impulse[i] - last_impulse[i]) / (time - last_time);

  // save the last condition
  last_time = time;
  for (size_t i=0; i<Dimensions; ++i) last_impulse[i] = this_impulse[i];

  return forces;
}

// Calculate total system impulse
std::array<float,Dimensions>
Simulation::calculate_total_impulse() {

  std::array<float,Dimensions> impulse = {0.0};

  // calculate impulse from particles
  for (auto &src : vort) {
    std::array<float,Dimensions> imp = std::visit([=](auto& elem) { return elem.get_total_impulse(); }, src);
    for (size_t i=0; i<Dimensions; ++i) impulse[i] += imp[i];
  }
  // then add up the impulse from bodies - DO WE NEED TO RE-SOLVE BEM FIRST?
  for (auto &src : bdry) {
    std::array<float,Dimensions> imp = std::visit([=](auto& elem) { return elem.get_total_impulse(); }, src);
    for (size_t i=0; i<Dimensions; ++i) impulse[i] += imp[i];
  }

  std::cout << "  total impulse " << impulse[0] << " " << impulse[1] << " " << impulse[2] << std::endl;

  return impulse;
}

// Add elements - any kind
void Simulation::add_elements(const ElementPacket<float> _elems,
                              const elem_t _et, const move_t _mt,
                              std::shared_ptr<Body> _bptr) {

  // skip out early if nothing's here
  if (_elems.nelem == 0) return;

  // now split on which Collection will receive this
  if (_et == active) {
    // it's active vorticity, add to vort
    file_elements(vort, _elems, active, _mt, _bptr);
    // in that routine, we will look for a match for move type, body pointer, and points/surfs/vols
  } else if (_et == reactive) {
    file_elements(bdry, _elems, reactive, _mt, _bptr);
  } else {
    file_elements(fldpt, _elems, inert, _mt, _bptr);
  }
}

// File the new elements into the correct collection
void Simulation::file_elements(std::vector<Collection>& _collvec,
                               const ElementPacket<float> _elems,
                               const elem_t _et, const move_t _mt,
                               std::shared_ptr<Body> _bptr) {

  // search the collections list for a match (same movement type, Body, elem dims)
  size_t imatch = 0;
  bool no_match = true;
  for (size_t i=0; i<_collvec.size(); ++i) {
    // assume match
    bool this_match = true;

    // check movement type
    const move_t tmt = std::visit([=](auto& elem) { return elem.get_movet(); }, _collvec[i]);
    if (_mt != tmt) {
      this_match = false;

    } else if (tmt == bodybound) {
      // check body pointer
      std::shared_ptr<Body> tbp = std::visit([=](auto& elem) { return elem.get_body_ptr(); }, _collvec[i]);
      if (_bptr != tbp) this_match = false;
    }

    // check Collections element dimension
    auto& coll = _collvec[i];
    if (std::holds_alternative<Points<float>>(coll) and _elems.ndim != 0) {
      this_match = false;
    } else if (std::holds_alternative<Surfaces<float>>(coll) and _elems.ndim != 2) {
      this_match = false;
    }

    if (this_match) {
      imatch = i;
      no_match = false;
    }
  }

  // if no match, or no collections exist
  if (no_match) {
    // make a new collection according to element dimension
    if (_elems.ndim == 0) {
      _collvec.push_back(Points<float>(_elems, _et, _mt, _bptr, get_vdelta()));
#ifdef USE_OGL_COMPUTE
      { // tell the new collection where the compute shader vao is
        Points<float>& pts = std::get<Points<float>>(_collvec.back());
        pts.set_opengl_compute_state(cgl);
      }
#endif
    } else if (_elems.ndim == 2) {
      _collvec.push_back(Surfaces<float>(_elems, _et, _mt, _bptr));
#ifdef USE_OGL_COMPUTE
      // tell it where the compute shader vao is
      {
        Surfaces<float>& surf = std::get<Surfaces<float>>(_collvec.back());
        surf.set_opengl_compute_state(cgl);
      }
#endif
    }

  } else {
    // found a match - get the Collection that matched
    auto& coll = _collvec[imatch];

    // proceed to add the correct object type
    if (_elems.ndim == 0) {
      Points<float>& pts = std::get<Points<float>>(coll);
      pts.add_new(_elems, get_vdelta());
    } else if (_elems.ndim == 2) {
      Surfaces<float>& surf = std::get<Surfaces<float>>(coll);
      surf.add_new(_elems);
    }
  }
}

// add a new Body with the given name
void Simulation::add_body(std::shared_ptr<Body> _body) {
  bodies.emplace_back(_body);
  std::cout << "  added new body (" << _body->get_name() << "), now have " << bodies.size() << std::endl;
}

// return a Body pointer to the last Body in the array
std::shared_ptr<Body> Simulation::get_last_body() {
  std::shared_ptr<Body> bp;

  if (bodies.size() == 0) {
    std::cout << "  no last body found, creating (ground)" << std::endl;
    bp = std::make_shared<Body>();
    bp->set_name("ground");
    add_body(bp);
  } else {
    bp = bodies.back();
    std::cout << "  returning last body (" << bp->get_name() << ")" << std::endl;
  }

  return bp;
}

// return a Body pointer to the body matching the given name
std::shared_ptr<Body> Simulation::get_pointer_to_body(const std::string _name) {
  std::shared_ptr<Body> bp;

  for (auto &bptr : bodies) {
    if (_name.compare(bptr->get_name()) == 0) {
      std::cout << "  found body matching name (" << _name << ")" << std::endl;
      bp = bptr;
    }
  }

  // or ground if none match
  if (not bp) {
    std::cout << "  no body matching (" << _name << ") found, creating (ground)" << std::endl;
    bp = std::make_shared<Body>();
    bp->set_name("ground");
    add_body(bp);
  }

  return bp;
}

// check vs. step and time to see if simulation should pause/stop
bool Simulation::test_vs_stop() {
  bool should_stop = false;
  if (using_max_steps() and get_max_steps() == nstep) {
    std::cout << "Stopping at step " << get_max_steps() << std::endl;
    should_stop = true;
  }
  if (using_end_time() and get_end_time() <= time+0.5*dt){
    std::cout << "Stopping at time " << get_end_time() << std::endl;
    should_stop = true;
  }
  return should_stop;
}

// this is different because we have to trigger when last step is still running
bool Simulation::test_vs_stop_async() {
  bool should_stop = false;
  static bool already_reported = false;

  if (using_max_steps() and get_max_steps() == nstep+1) {
    if (not already_reported) {
      std::cout << std::endl << "Stopping at step " << get_max_steps() << std::endl;
      already_reported = true;
    }
    should_stop = true;
  }

  if (using_end_time() and get_end_time() >= time+0.5*dt
                       and get_end_time() <= time+1.5*dt) {
    if (not already_reported) {
      std::cout << std::endl << "Stopping at time " << get_end_time() << std::endl;
      already_reported = true;
    }
    should_stop = true;
  }

  // reset the toggle
  if (not should_stop) already_reported = false;

  return should_stop;
}

