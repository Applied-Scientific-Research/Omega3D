/*
 * main_gui.cpp - Driver code for Omega3D + ImGui + Vc vortex particle method
 *                and boundary element method solver, GUI version
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "Simulation.h"
#include "JsonHelper.h"
#include "Body.h"
#include "RenderParams.h"

#ifdef _WIN32
  // for glad
  #define APIENTRY __stdcall
  // for C++11 stuff
  #include <ciso646>
#endif
#include "glad.h"

// header-only immediate-mode GUI
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw_gl3.h"
#include "imgui/ImguiWindowsFileIO.hpp"

// header-only png writing
#include "stb/FrameBufferToImage.h"

//#include <GL/gl3w.h>    // This example is using gl3w to access OpenGL
// functions (because it is small). You may use glew/glad/glLoadGen/etc.
// whatever already works for you.
#include <GLFW/glfw3.h>

#include <cstdio>
#include <iostream>
#include <vector>

static void error_callback(int error, const char* description) {
  fprintf(stderr, "Error %d: %s\n", error, description);
}

//static void keyboard_callback(int key, int action) {
//  printf("%d %d\n", key, action);
//}

//
// this is NOT a GLFW callback, but my own, and it needs the
// ImGuiIO data structure for information on the mouse state
//
void mouse_callback(GLFWwindow* /*_thiswin*/,
                    ImGuiIO& io,
                    float*   _cx,
                    float*   _cy,
                    float*   _size) {

  // first, use left-click and drag to move the data
  static bool lbutton_down = false;

  if (io.MouseClicked[0]) lbutton_down = true;
  if (io.MouseReleased[0]) lbutton_down = false;

  if (lbutton_down) {
    //std::cout << "free mouse moved " << io.MouseDelta.x << " " << io.MouseDelta.y << std::endl;
    // do your drag here

    // this worked on Linux:
    //int display_w, display_h;
    //glfwGetFramebufferSize(_thiswin, &display_w, &display_h);
    //(*_cx) -= 2.0f * (*_size) * (float)io.MouseDelta.x / (float)display_w;
    //(*_cy) += 2.0f * (*_size) * (float)io.MouseDelta.y / (float)display_w;

    // this works on a Retina display:
    (*_cx) -= 2.0f * (*_size) * (float)io.MouseDelta.x / io.DisplaySize.x;
    (*_cy) += 2.0f * (*_size) * (float)io.MouseDelta.y / io.DisplaySize.x;
  }

  // then, use scroll wheel to zoom!
  //std::cout << "free mouse wheel " << io.MouseWheel << std::endl;
  if (io.MouseWheel != 0) {
    // do your drag here
    //int display_w, display_h;
    //glfwGetFramebufferSize(_thiswin, &display_w, &display_h);

    // change the size
    const float oldsize = (*_size);
    (*_size) *= std::pow(1.1f, io.MouseWheel);

    // and adjust the center such that the zoom occurs about the pointer!
    //const float ar = (float)display_h / (float)display_w;
    const float ar = io.DisplaySize.y / io.DisplaySize.x;

    // this only scales around world origin
    //(*_cx) += 2.0f * ((float)io.MousePos.x / (float)display_w - 0.5f) * (oldsize - (*_size));
    //(*_cy) += 2.0f * (0.5f - (float)io.MousePos.y / (float)display_h) * (oldsize - (*_size)) * ar;
    (*_cx) += 2.0f * ((float)io.MousePos.x / io.DisplaySize.x - 0.5f) * (oldsize - (*_size));
    (*_cy) += 2.0f * (0.5f - (float)io.MousePos.y / io.DisplaySize.y) * (oldsize - (*_size)) * ar;
  }
}

//
// Helper routine to determine orthographic projection matrix
// given coords at screen center and a measure of size
// Also changes overall pixels-to-length scale
//
void compute_ortho_proj_mat(GLFWwindow*         _thiswin,
                            const float         _cx,
                            const float         _cy,
                            float*              _size,
                            std::vector<float>& _projmat) {

  // track changes in window!
  static int last_w, last_h = -1;

  // get current window size
  int display_w, display_h;
  glfwGetFramebufferSize(_thiswin, &display_w, &display_h);

  // compare window size to previous call
  if (last_h != -1) {
    if (last_h != display_h or last_w != display_w) {
      // window aspect ratio changed, adjust _size
      (*_size) *= sqrt(  ((float)last_h   /(float)last_w   )
                       / ((float)display_h/(float)display_w));
    }
  }

  const float vsx = (*_size);
  const float vsy = (*_size) * (float)display_h / (float)display_w;
  _projmat =
    { 1.0f/vsx, 0.0f,     0.0f, 0.0f,
      0.0f,     1.0f/vsy, 0.0f, 0.0f,
      0.0f,     0.0f,    -1.0f, 0.0f,
     -_cx/vsx, -_cy/vsy,  0.0f, 1.0f };

  // save window size for next call
  last_w = display_w;
  last_h = display_h;
}

//
// resize a window and framebuffer programmatically
//
void resize_to_resolution(GLFWwindow* window, const int new_w, const int new_h) {

  // get framebuffer size
  int fb_w, fb_h;
  glfwGetFramebufferSize(window, &fb_w, &fb_h);
  //std::cout << "Framebuffer size is " << fb_w << " x " << fb_h << std::endl;

  // get window size
  int ws_w, ws_h;
  glfwGetWindowSize(window, &ws_w, &ws_h);
  //std::cout << "Window size is " << ws_w << " x " << ws_h << std::endl;

  // on normal monitors, these numbers should be the same; on retina displays, they may not

  // check and resize
  if (fb_w != new_w or fb_h != new_h) {
    // on a retina display...do anything different?

    glfwSetWindowSize(window, new_w, new_h);
    std::cout << "Resizing window/framebuffer to " << new_w << " x " << new_h << std::endl;
  }
}

static void ShowHelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(450.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}


// execution starts here

int main(int argc, char const *argv[]) {
  std::cout << std::endl << "Omega3D GUI" << std::endl;

  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;
  std::vector< std::unique_ptr<MeasureFeature> > mfeatures;
  size_t nframes = 0;
  static bool sim_is_running = false;
  static bool begin_single_step = false;

  // Set up primary OpenGL window
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    return 1;
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
  GLFWwindow* window = glfwCreateWindow(1280, 720, "Omega3D GUI", nullptr, nullptr);
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  //gl3wInit();

  if (!gladLoadGL()) {
    printf("Something went wrong!\n");
    exit(-1);
  }

  // Setup ImGui binding
  ImGui_ImplGlfwGL3_Init(window, true);

  //glfwSetKeyCallback(keyboard_callback);

  //glfwSetWindowCloseCallback(window, window_close_callback);

  // Load Fonts
  // (there is a default font, this is only if you want to change it. see extra_fonts/README.txt for more details)

  // Get and set some IO functions
  ImGuiIO& io = ImGui::GetIO();
  io.IniFilename = ".omega3d.ini";
  std::vector<std::string> recent_json_files;
  std::vector<std::string> recent_geom_files;

  // a string to hold any error messages
  std::string sim_err_msg;

  // GUI and drawing parameters
  bool export_vtk_this_frame = false;	// write a vtk with the current data
  bool draw_this_frame = false;		// draw the frame immediately
  bool record_all_frames = false;	// save a frame when a new one is ready
  bool show_stats_window = false;
  bool show_terminal_window = false;
  bool show_test_window = false;
  bool show_geom_input_window = false;
  bool show_json_input_window = false;
  bool show_file_output_window = false;
  bool show_bdry_create_window = false;
  //static bool show_origin = true;
  static bool is_viscous = false;

  // colors and projection matrix for the render view
  RenderParams rparams;
  std::vector<float> gl_projection;
  compute_ortho_proj_mat(window, rparams.vcx, rparams.vcy, &rparams.vsize, gl_projection);


  // Main loop
  while (!glfwWindowShouldClose(window))
  {
    glfwPollEvents();
    ImGui_ImplGlfwGL3_NewFrame();

    //
    // Initialize simulation
    //

    if (not sim.is_initialized() and (sim_is_running || begin_single_step)) {

      std::cout << std::endl << "Initializing simulation" << std::endl;

      // initialize particle distributions
      for (auto const& ff: ffeatures) {
        sim.add_particles( ff->init_particles(sim.get_ips()) );
      }

      // initialize solid objects
      for (auto const& bf : bfeatures) {
        sim.add_boundary( bf->get_body(), bf->init_elements(sim.get_ips()) );
      }

      // initialize measurement features
      for (auto const& mf: mfeatures) {
        sim.add_fldpts( mf->init_particles(rparams.tracer_scale*sim.get_ips()), mf->moves() );
      }

      sim.set_initialized();

      // check setup for obvious errors
      sim_err_msg = sim.check_initialization();

      if (sim_err_msg.empty()) {
        // begin the initial calculation at t=0 so we can draw something
        sim.async_first_step();

      } else {
        // the last step had some difficulty
        std::cout << std::endl << "ERROR: " << sim_err_msg;
        // stop the run
        sim_is_running = false;
      }

      begin_single_step = false;
    }

    //
    // Update simulation
    //

    // get results of latest step, if it just completed
    bool is_ready = sim.test_for_new_results();

    // before we start again, write the vtu output
    if (sim.get_nparts() > 0 and export_vtk_this_frame) {
      sim.write_vtk();
      export_vtk_this_frame = false;
    }

    // see if we should start a new step
    if (is_ready and (sim_is_running || begin_single_step)) {

      // check flow for blow-up or dynamic errors
      sim_err_msg = sim.check_simulation();

      if (sim_err_msg.empty()) {
        // the last simulation step was fine, OK to continue

        // generate new particles from emitters
        for (auto const& ff : ffeatures) {
          sim.add_particles( ff->step_particles(sim.get_ips()) );
        }
        for (auto const& mf: mfeatures) {
          sim.add_fldpts( mf->step_particles(rparams.tracer_scale*sim.get_ips()), mf->moves() );
        }

        // begin a new dynamic step: convection and diffusion
        sim.async_step();

      } else {
        // the last step had some difficulty
        std::cout << std::endl << "ERROR: " << sim_err_msg;

        // stop the run
        sim_is_running = false;
      }

      begin_single_step = false;
    }

    // check the error message and display if there is one
    if (not sim_err_msg.empty()) {
      // write a warning/error message
      ImGui::OpenPopup("Simulation error occurred");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("Simulation error occurred")) {
        ImGui::Spacing();
        ImGui::TextWrapped(sim_err_msg.c_str());
        ImGui::Spacing();
        if (ImGui::Button("Got it.", ImVec2(120,0))) {
          // clear out the error message first
          sim_err_msg.clear();
          ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
      }
    }


    // check mouse for drag and rescaling!
    if (not io.WantCaptureMouse) {
      mouse_callback(window, io, &rparams.vcx, &rparams.vcy, &rparams.vsize);
    }

    // check for keypresses to toggle state
    //if (not io.WantCaptureKeyboard) {
      //keyboard_callback(
    //}

    //
    // The main Omega3D window
    //
    {

    ImGui::Begin("Omega3D");
    ImGui::TextWrapped("Welcome to Omega3D. Select a simulation from the drop-down, or manually set your simulation globals, then add one or more flow structures. Space bar starts and stops the run. Have fun!");
    ImGui::Spacing();

    // Select pre-populated simulations
    {
      static int sim_item = 0;
      const char* sim_items[] = { "Select a simulation...", "single vortex ring - no viscosity", "leapfrogging vortex rings - no viscosity", "colliding vortex rings - no viscosity", "single viscous vortex ring" };
      ImGui::Combo("", &sim_item, sim_items, 5);

      float* dt = sim.addr_dt();
      float* fs = sim.addr_fs();
      float* re = sim.addr_re();
      std::shared_ptr<Body> bp;

      switch(sim_item) {
        case 0:
          // nothing
          break;
        case 1:
          // one singular vortex ring
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.002;
          fs[0] = 0.0; fs[1] = 0.0;
          sim.set_re_for_ips(0.015);
          // generate the features
          ffeatures.emplace_back(std::make_unique<SingularRing>(0.0,0.0,0.0, 0.9,0.05,0.1, 0.5, 1.0));
          is_viscous = false;
          sim.set_diffuse(false);
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
        case 2:
          // leapfrogging vortex rings
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.002;
          fs[0] = 0.0; fs[1] = 0.0;
          sim.set_re_for_ips(0.02);
          // generate the features
          ffeatures.emplace_back(std::make_unique<SingularRing>(0.0,0.0,0.0, 0.9,0.05,0.1, 0.5, 1.0));
          ffeatures.emplace_back(std::make_unique<SingularRing>(-0.18,-0.01,-0.02, 0.9,0.05,0.1, 0.5, 1.0));
          is_viscous = false;
          sim.set_diffuse(false);
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
        case 3:
          // colliding vortex rings
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.002;
          fs[0] = 0.0; fs[1] = 0.0;
          sim.set_re_for_ips(0.02);
          // generate the features
          ffeatures.emplace_back(std::make_unique<SingularRing>(0.4,0.0,0.0, -0.9,-0.05,-0.1, 0.5, 1.0));
          ffeatures.emplace_back(std::make_unique<SingularRing>(0.04,-0.02,-0.04, 0.9,0.05,0.1, 0.5, 1.0));
          is_viscous = false;
          sim.set_diffuse(false);
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
        case 4:
          // Re=100 vortex ring
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.01;
          fs[0] = 0.0; fs[1] = 0.0;
          *re = 100.0;
          // generate the features
          //ffeatures.emplace_back(std::make_unique<SingularRing>(0.0,0.0,0.0, 1.0,0.0,0.0, 0.5, 1.0));
          ffeatures.emplace_back(std::make_unique<SingularRing>(0.0,0.0,0.0, 0.9,0.05,0.1, 0.5, 1.0));
          is_viscous = true;
          sim.set_diffuse(true);
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
      } // end switch
    }

    // or load a simulation from a JSON file
    ImGui::SameLine();
    if (ImGui::Button("Or load a json file", ImVec2(160,0))) show_json_input_window = true;

    if (show_json_input_window) {
      bool try_it = false;
      static std::string infile = "input.json";

      if (fileIOWindow( try_it, infile, recent_json_files, "Open", {"*.json", "*.*"}, true, ImVec2(500,250))) {
        show_json_input_window = false;

        if (try_it and !infile.empty()) {
          // remember
          recent_json_files.push_back( infile );

          // stop and clear before loading
          sim.reset();
          bfeatures.clear();
          ffeatures.clear();

          // load and report
          read_json(sim, ffeatures, bfeatures, mfeatures, rparams, infile);

          // we have to manually set this variable
          is_viscous = sim.get_diffuse();
          // run one step so we know what we have
          begin_single_step = true;

          // check and possibly resize the window to match the saved resolution
          resize_to_resolution(window, rparams.width, rparams.height);
        }
      }
    }
    ImGui::Spacing();


    //if (ImGui::CollapsingHeader("Simulation globals", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (ImGui::CollapsingHeader("Simulation globals")) {

      ImGui::SliderFloat("Time step", sim.addr_dt(), 0.001f, 0.1f, "%.4f", 2.0f);
      ImGui::Checkbox("Fluid is viscous (diffuses)", &is_viscous);
      if (is_viscous) {
        // show the toggle for AMR
        //static bool use_amr = false;
        //ImGui::Checkbox("Allow adaptive resolution", &use_amr);
        //sim.set_amr(use_amr);
        sim.set_diffuse(true);
        // and let user choose Reynolds number
        ImGui::SliderFloat("Reynolds number", sim.addr_re(), 10.0f, 2000.0f, "%.1f", 2.0f);
        ImGui::Text("Particle spacing %g", sim.get_ips());
      } else {
        static float my_ips = 0.03141;
        ImGui::SliderFloat("Particle spacing", &my_ips, 0.001f, 0.1f, "%.3f", 2.0f);
        // change underlying Re when this changes
        sim.set_re_for_ips(my_ips);
        my_ips = sim.get_ips();
      }
      ImGui::InputFloat3("Freestream speed", sim.addr_fs());

      // set stop/pause conditions
      bool use_step_pause = sim.using_max_steps();
      ImGui::Checkbox("Pause at step", &use_step_pause);
      if (use_step_pause) {
        int pause_step = sim.get_max_steps();
        ImGui::SameLine();
        ImGui::InputInt(" ", &pause_step);
        sim.set_max_steps((size_t)pause_step);
      } else {
        sim.unset_max_steps();
      }
      bool use_time_pause = sim.using_end_time();
      ImGui::Checkbox("Pause at time", &use_time_pause);
      if (use_time_pause) {
        float pause_time = sim.get_end_time();
        ImGui::SameLine();
        ImGui::InputFloat(" ", &pause_time);
        //ImGui::InputFloat("Pause at time", float* v, float step = 0.0f, float step_fast = 0.0f, int decimal_precision = -1, ImGuiInputTextFlags extra_flags = 0);
        sim.set_end_time((double)pause_time);
      } else {
        sim.unset_end_time();
      }
    }

    ImGui::Spacing();
    //if (ImGui::CollapsingHeader("Flow structures", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (ImGui::CollapsingHeader("Flow structures")) {

      int buttonIDs = 10;

      // list existing flow features here
      int del_this_item = -1;
      for (int i=0; i<(int)ffeatures.size(); ++i) {
        // add a "remove" button here somehow
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_item = i;
        ImGui::PopID();

        //ImGui::SameLine();
        //ImGui::PushID(++buttonIDs);
        //if (ImGui::SmallButton("edit", ImVec2(60,0))) del_this_item = 0;
        //ImGui::PopID();

        ImGui::SameLine();
        ImGui::Text("%s", ffeatures[i]->to_string().c_str());
      }
      if (del_this_item > -1) {
        std::cout << "Asked to delete flow feature " << del_this_item << std::endl;
        ffeatures.erase(ffeatures.begin()+del_this_item);
      }

      // list existing boundary features here
      int del_this_bdry = -1;
      for (int i=0; i<(int)bfeatures.size(); ++i) {
        // add a "remove" button here somehow
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_bdry = i;
        ImGui::PopID();

        ImGui::SameLine();
        ImGui::Text("%s", bfeatures[i]->to_string().c_str());
      }
      if (del_this_bdry > -1) {
        std::cout << "Asked to delete boundary feature " << del_this_bdry << std::endl;
        bfeatures.erase(bfeatures.begin()+del_this_bdry);
      }

      // list existing measurement features here
      int del_this_measure = -1;
      for (int i=0; i<(int)mfeatures.size(); ++i) {
        // add a "remove" button here somehow
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_measure = i;
        ImGui::PopID();

        ImGui::SameLine();
        ImGui::Text("%s", mfeatures[i]->to_string().c_str());
      }
      if (del_this_measure > -1) {
        std::cout << "Asked to delete measurement feature " << del_this_measure << std::endl;
        mfeatures.erase(mfeatures.begin()+del_this_measure);
      }

      // button and modal window for adding new flow structures
      if (ImGui::Button("Add flow structure")) ImGui::OpenPopup("New flow structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New flow structure"))
      {
        static int item = 1;
        const char* items[] = { "vortex blob", "random particles", "singular vortex ring", "thick vortex ring" };
        ImGui::Combo("type", &item, items, 4);

        static float xc[3] = {0.0f, 0.0f, 0.0f};	// a center
        //static float rad = 5.0 * sim.get_ips();	// a major radius
        static float rad = 1.0;				// a major radius
        static float soft = sim.get_ips();		// a softness or minor radius
        static float vstr[3] = {0.0f, 0.0f, 1.0f};	// a vectorial strength
        static float strmag = 1.0f;			// a scalar strength
        static float circ = 1.0f;			// a circulation
        static int npart = 1000;
        static float xs[3] = {2.0f, 2.0f, 2.0f};	// a size
        int guess_n = 0;

        // always ask for center
        ImGui::InputFloat3("center", xc);

        // show different inputs based on what is selected
        switch(item) {
          case 0:
            // a blob of multiple vortons
            ImGui::InputFloat3("strength", vstr);
            ImGui::SliderFloat("radius", &rad, sim.get_ips(), 10.0f*sim.get_ips(), "%.4f");
            ImGui::SliderFloat("softness", &soft, sim.get_ips(), rad, "%.4f");
            guess_n = 4.1888f * std::pow((2.0f*rad+soft)/sim.get_ips(), 3);
            ImGui::TextWrapped("This feature will add about %d particles", guess_n);
            if (ImGui::Button("Add vortex blob")) {
              ffeatures.emplace_back(std::make_unique<VortexBlob>(xc[0],xc[1],xc[2], vstr[0],vstr[1],vstr[2], rad, soft));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            // it would be nice to be able to put this all in
            //SingleParticle::draw_creation_gui();
            break;

          case 1:
            // random particles in a block
            ImGui::SliderInt("number", &npart, 10, 100000);
            ImGui::SliderFloat3("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
            ImGui::SliderFloat("strength magnitude", &strmag, 0.01f, 10.0f, "%.3f", 2.0f);
            ImGui::TextWrapped("This feature will add %d particles", npart);
            if (ImGui::Button("Add random vorticies")) {
              ffeatures.emplace_back(std::make_unique<BlockOfRandom>(xc[0],xc[1],xc[2], xs[0],xs[1],xs[2], strmag, npart));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;

          case 2:
            // oriented singular vortex ring
            ImGui::InputFloat3("direction", vstr);
            ImGui::SliderFloat("circulation", &circ, 0.001f, 10.0f, "%.3f");
            ImGui::SliderFloat("radius", &rad, 3.0f*sim.get_ips(), 10.0f, "%.3f");
            guess_n = 1 + (2.0f * 3.1416f * rad / sim.get_ips());
            ImGui::TextWrapped("This feature will add about %d particles", guess_n);
            if (ImGui::Button("Add singular vortex ring")) {
              ffeatures.emplace_back(std::make_unique<SingularRing>(xc[0],xc[1],xc[2], vstr[0],vstr[1],vstr[2], rad, circ));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;

          case 3:
            // oriented thick-cored vortex ring
            ImGui::InputFloat3("direction", vstr);
            ImGui::SliderFloat("circulation", &circ, 0.001f, 10.0f, "%.4f");
            ImGui::SliderFloat("radius", &rad, 3.0f*sim.get_ips(), 10.0f, "%.3f");
            ImGui::SliderFloat("thickness", &soft, sim.get_ips(), 10.0f*sim.get_ips(), "%.4f");
            guess_n = (1 + (2.0f * 3.1416f * rad / sim.get_ips())) * std::pow(soft / sim.get_ips(), 2);
            ImGui::TextWrapped("This feature will add about %d particles", guess_n);
            if (ImGui::Button("Add thick vortex ring")) {
              ffeatures.emplace_back(std::make_unique<ThickRing>(xc[0],xc[1],xc[2], vstr[0],vstr[1],vstr[2], rad, soft, circ));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      }


      // button and modal window for adding new boundary objects
      ImGui::SameLine();
      if (ImGui::Button("Add boundary structure")) show_bdry_create_window = true;
      if (show_bdry_create_window)
      {
        ImGui::SetNextWindowSize(ImVec2(400,275), ImGuiSetCond_FirstUseEver);
        ImGui::Begin("New boundary structure", &show_bdry_create_window);

        // define movement first
        static int mitem = 0;
        const char* mitems[] = { "fixed", "attached to previous", "according to formula" };
        //const char* mitems[] = { "fixed", "attached to previous", "according to formula", "dynamic" };
        ImGui::Combo("movement", &mitem, mitems, 3);
        static char strx[512] = "0.0*t";
        static char stry[512] = "0.0*t";
        static char strz[512] = "0.0*t";

        // show different inputs based on what is selected
        switch(mitem) {
          case 0:
            // this geometry is fixed (attached to inertial)
            break;
          case 1:
            // this geometry is attached to the previous geometry
            break;
          case 2:
            // this geometry is attached to a new moving body
            ImGui::InputText("x position", strx, 512);
            ImGui::SameLine();
            ShowHelpMarker("Use C-style expressions, t is time\n+ - / * \% ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
            ImGui::InputText("y position", stry, 512);
            ImGui::SameLine();
            ShowHelpMarker("Use C-style expressions, t is time\n+ - / * \% ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
            ImGui::InputText("z position", stry, 512);
            ImGui::SameLine();
            ShowHelpMarker("Use C-style expressions, t is time\n+ - / * \% ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
            break;
        }

        // define geometry second
        static int item = 0;
        const char* items[] = { "sphere", "rectangle", "from file" };
        ImGui::Combo("type", &item, items, 3);

        static float xc[3] = {0.0f, 0.0f, 0.0f};
        //static float rotdeg = 0.0f;
        static float scale = 1.0;

        // always ask for center
        ImGui::InputFloat3("center", xc);

        // show different inputs based on what is selected
        switch(item) {
          case 0:
            // create a solid sphere
            ImGui::SliderFloat("diameter", &scale, 0.01f, 10.0f, "%.4f", 2.0);
            ImGui::TextWrapped("This feature will add a solid spherical body centered at the given coordinates");
            if (ImGui::Button("Add spherical body")) {
              std::shared_ptr<Body> bp;
              switch(mitem) {
                case 0:
                  // this geometry is fixed (attached to inertial)
                  bp = sim.get_pointer_to_body("ground");
                  break;
                case 1:
                  // this geometry is attached to the previous geometry (or ground)
                  bp = sim.get_last_body();
                  break;
                case 2:
                  // this geometry is attached to a new moving body
                  bp = std::make_shared<Body>();
                  bp->set_pos(0, std::string(strx));
                  bp->set_pos(1, std::string(stry));
                  bp->set_pos(2, std::string(strz));
                  bp->set_name("sphere");
                  sim.add_body(bp);
                  break;
              }
              //bfeatures.emplace_back(std::make_unique<SolidCircle>(bp, xc[0], xc[1], xc[2], scale));
              //std::cout << "Added " << (*bfeatures.back()) << std::endl;
              show_bdry_create_window = false;
            }
            ImGui::SameLine();
            break;
          case 1:
            // create a solid rectangle
            static float xs[3] = {1.0f, 1.0f, 1.0f};
            ImGui::InputFloat3("side lengths", xs);
            //ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 89.0f, "%.0f");
            //ImGui::SliderAngle("orientation", &rotdeg);
            ImGui::TextWrapped("This feature will add a solid rectangular body centered at the given coordinates");
            if (ImGui::Button("Add rectangular body")) {
              std::shared_ptr<Body> bp;
              switch(mitem) {
                case 0:
                  // this geometry is fixed (attached to inertial)
                  bp = sim.get_pointer_to_body("ground");
                  break;
                case 1:
                  // this geometry is attached to the previous geometry (or ground)
                  bp = sim.get_last_body();
                  break;
                case 2:
                  // this geometry is attached to a new moving body
                  bp = std::make_shared<Body>();
                  bp->set_pos(0, std::string(strx));
                  bp->set_pos(1, std::string(stry));
                  bp->set_pos(2, std::string(strz));
                  bp->set_name("rectangle");
                  sim.add_body(bp);
                  break;
              }
              //bfeatures.emplace_back(std::make_unique<SolidSquare>(bp, xc[0], xc[1], sqside, rotdeg));
              //std::cout << "Added " << (*bfeatures.back()) << std::endl;
              show_bdry_create_window = false;
            }
            ImGui::SameLine();
            break;
          case 2:
            // load a geometry file
            static std::string infile = "input.obj";
            static std::string shortname = infile;
            if (ImGui::Button("Geometry file", ImVec2(200,0))) show_geom_input_window = true;

            if (show_geom_input_window) {
              bool try_it = false;

              if (fileIOWindow( try_it, infile, recent_geom_files, "Open", {"*.obj", "*.stl", "*.ply", "*.*"}, true, ImVec2(500,250))) {
                show_geom_input_window = false;

                if (try_it and !infile.empty()) {
                  // remember
                  recent_geom_files.push_back( infile );

                  // now remove the leading directories from the string
                  const size_t lastchar = infile.find_last_of("/\\");
                  shortname = infile.substr(lastchar+1);
                }
              }
            }
            ImGui::SameLine();
            ImGui::Text(shortname.c_str());
            ImGui::SliderFloat("diameter", &scale, 0.01f, 10.0f, "%.4f", 2.0);
            ImGui::TextWrapped("This feature will add a solid body centered at the given coordinates");
            if (ImGui::Button("Add geometry from file")) {
              std::shared_ptr<Body> bp;
              switch(mitem) {
                case 0:
                  // this geometry is fixed (attached to inertial)
                  bp = sim.get_pointer_to_body("ground");
                  break;
                case 1:
                  // this geometry is attached to the previous geometry (or ground)
                  bp = sim.get_last_body();
                  break;
                case 2:
                  // this geometry is attached to a new moving body
                  bp = std::make_shared<Body>();
                  bp->set_pos(0, std::string(strx));
                  bp->set_pos(1, std::string(stry));
                  bp->set_pos(2, std::string(strz));
                  bp->set_name("geom from file");
                  sim.add_body(bp);
                  break;
              }
              bfeatures.emplace_back(std::make_unique<ExteriorFromFile>(bp, xc[0], xc[1], xc[2], scale, scale, scale, infile));
              std::cout << "Added " << (*bfeatures.back()) << std::endl;
              show_bdry_create_window = false;
            }
            ImGui::SameLine();
            break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { show_bdry_create_window = false; }
        ImGui::End();
      }


      // button and modal window for adding new measurement objects
      ImGui::SameLine();
      if (ImGui::Button("Add measurement structure")) ImGui::OpenPopup("New measurement structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New measurement structure"))
      {
        static int item = 0;
        const char* items[] = { "single point/tracer", "streakline", "circle of tracers", "line of tracers", "measurement line" };
        ImGui::Combo("type", &item, items, 5);

        static float xc[3] = {0.0f, 0.0f, 0.0f};
        static float xf[3] = {0.0f, 1.0f, 0.0f};
        static bool is_lagrangian = true;
        static float rad = 2.0 * sim.get_ips();

        // show different inputs based on what is selected
        switch(item) {
          case 0:
            // a single measurement point
            ImGui::InputFloat3("position", xc);
            ImGui::Checkbox("Point follows flow", &is_lagrangian);
            ImGui::TextWrapped("This feature will add 1 point");
            if (ImGui::Button("Add single point")) {
              mfeatures.emplace_back(std::make_unique<SinglePoint>(xc[0], xc[1], xc[2], is_lagrangian));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
          case 1:
            // a tracer emitter
            ImGui::InputFloat3("position", xc);
            ImGui::TextWrapped("This feature will add 1 tracer emitter");
            if (ImGui::Button("Add streakline")) {
              mfeatures.emplace_back(std::make_unique<TracerEmitter>(xc[0], xc[1], xc[2]));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
          case 2:
            // a tracer circle
            ImGui::InputFloat3("center", xc);
            ImGui::SliderFloat("radius", &rad, 0.5f*sim.get_ips(), 0.5f, "%.4f");
            ImGui::TextWrapped("This feature will add about %d field points",
                               (int)(0.6*std::pow(2*rad/(rparams.tracer_scale*sim.get_ips()), 3)));
            if (ImGui::Button("Add circle of tracers")) {
              mfeatures.emplace_back(std::make_unique<TracerBlob>(xc[0], xc[1], xc[2], rad));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
          case 3:
            // a tracer line
            ImGui::InputFloat3("start", xc);
            ImGui::InputFloat3("finish", xf);
            ImGui::TextWrapped("This feature will add about %d field points",
                               1+(int)(std::sqrt(std::pow(xf[0]-xc[0],2)+std::pow(xf[1]-xc[1],2)+std::pow(xf[2]-xc[2],2))/(rparams.tracer_scale*sim.get_ips())));
            if (ImGui::Button("Add line of tracers")) {
              mfeatures.emplace_back(std::make_unique<TracerLine>(xc[0], xc[1], xc[2], xf[0], xf[1], xf[2]));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
          case 4:
            // a static, measurement line
            ImGui::InputFloat3("start", xc);
            ImGui::InputFloat3("finish", xf);
            ImGui::TextWrapped("This feature will add about %d field points",
                               1+(int)(std::sqrt(std::pow(xf[0]-xc[0],2)+std::pow(xf[1]-xc[1],2)+std::pow(xf[2]-xc[2],2))/(rparams.tracer_scale*sim.get_ips())));
            if (ImGui::Button("Add line of measurement points")) {
              mfeatures.emplace_back(std::make_unique<MeasurementLine>(xc[0], xc[1], xc[2], xf[0], xf[1], xf[2]));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      }

    } // end structure entry

    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Rendering parameters")) {
      ImGui::ColorEdit3("positive circulation", rparams.pos_circ_color);
      ImGui::ColorEdit3("negative circulation", rparams.neg_circ_color);
      ImGui::ColorEdit3("feature color",        rparams.default_color);
      ImGui::ColorEdit3("background color",     rparams.clear_color);
      //ImGui::Checkbox("show origin", &show_origin);
      ImGui::SliderFloat("vorticity density", &(rparams.circ_density), 0.01f, 2.0f, "%.2f", 2.0f);

      if (ImGui::Button("Recenter")) {
        // put everything back to center
        rparams.vcx = -0.5f;
        rparams.vcy = 0.0f;
        rparams.vsize = 2.0f;
      }

      // add button to recenter on all vorticity?

      // help separate this from the final info
      ImGui::Separator();
    }
    ImGui::Spacing();

    nframes++;

    // check vs. end conditions, if present
    if (sim.test_vs_stop_async()) sim_is_running = false;

    // all the other stuff
    {
      if (sim_is_running) {
        ImGui::Text("Simulation is running...step = %ld, time = %g", sim.get_nstep(), sim.get_time());
        if (ImGui::Button("PAUSE", ImVec2(200,0))) sim_is_running = false;
        // space bar pauses
        if (ImGui::IsKeyPressed(32)) sim_is_running = false;
      } else {
        ImGui::Text("Simulation is not running, step = %ld, time = %g", sim.get_nstep(), sim.get_time());
        if (ImGui::Button("RUN", ImVec2(200,0))) sim_is_running = true;
        ImGui::SameLine();
        if (ImGui::Button("Step", ImVec2(120,0))) begin_single_step = true;
        // and space bar resumes
        if (ImGui::IsKeyPressed(32)) sim_is_running = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Reset", ImVec2(120,0))) {
        std::cout << std::endl << "Reset requested" << std::endl;
        // remove all particles and reset timer
        sim.reset();
        std::cout << "Reset complete" << std::endl;
      }

      ImGui::Spacing();
      //ImGui::Separator();
      /*
      ImGui::Text("Open additional windows");
      if (ImGui::Button("Plot statistics")) show_stats_window ^= 1;
      ImGui::SameLine();
      if (ImGui::Button("Show terminal output")) show_terminal_window ^= 1;
      ImGui::SameLine();
      */
      //if (ImGui::Button("ImGui Samples")) show_test_window ^= 1;

      ImGui::Text("Draw frame rate: %.2f ms/frame (%.1f FPS)",
                  1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

      ImGui::Text("Number of panels: %ld  particles: %ld  field points: %ld",
                  sim.get_npanels(), sim.get_nparts(), sim.get_nfldpts());

      // save the simulation to a JSON or VTK file
      ImGui::Spacing();
      if (ImGui::Button("Save setup to json", ImVec2(180,0))) show_file_output_window = true;
      ImGui::SameLine();
      if (ImGui::Button("Save parts to vtk", ImVec2(170,0))) export_vtk_this_frame = true;

      if (show_file_output_window) {
        bool try_it = false;
        static std::string outfile = "output.json";

        if (fileIOWindow( try_it, outfile, recent_json_files, "Save", {"*.json", "*.*"}, false, ImVec2(500,250))) {
          show_file_output_window = false;

          if (try_it) {
            // remember
            recent_json_files.push_back( outfile );

            // retrieve window sizes
            glfwGetWindowSize(window, &rparams.width, &rparams.height);

            // write and echo
            write_json(sim, ffeatures, bfeatures, mfeatures, rparams, outfile);
            std::cout << std::endl << "Wrote simulation to " << outfile << std::endl;
          }
        }
      }

      // PNG output of the render frame
      if (ImGui::Button("Save screenshot to png", ImVec2(180,0))) draw_this_frame = true;
      ImGui::SameLine();
      if (record_all_frames) {
        if (ImGui::Button("STOP", ImVec2(80,0))) {
          record_all_frames = false;
          sim_is_running = false;
        }
      } else {
        if (ImGui::Button("RECORD to png", ImVec2(120,0))) {
          record_all_frames = true;
          sim_is_running = true;
        }
      }
    }

    // done drawing the UI window
    ImGui::End();
    }

    // Show the simulation stats as 2D plots
    if (show_stats_window)
    {
      ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiSetCond_FirstUseEver);
      ImGui::Begin("Statistics", &show_stats_window);
      ImGui::Text("Hello");
      ImGui::End();
    }

    // Show the terminal output of the program
    if (show_terminal_window)
    {
      ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiSetCond_FirstUseEver);
      ImGui::Begin("Terminal", &show_terminal_window);
      ImGui::Text("Hello");
      ImGui::End();
    }

    // Show the ImGui test window. Most of the sample code is in ImGui::ShowTestWindow()
    if (show_test_window)
    {
      ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
      ImGui::ShowTestWindow(&show_test_window);
    }

    // Rendering
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(rparams.clear_color[0], rparams.clear_color[1], rparams.clear_color[2], rparams.clear_color[3]);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw the simulation: panels and particles
    compute_ortho_proj_mat(window, rparams.vcx, rparams.vcy, &rparams.vsize, gl_projection);
    sim.drawGL(gl_projection, rparams);

    // if simulation has not been initted, draw the features instead!
    //for (auto const& bf : bfeatures) { bf.drawGL(gl_projection, rparams); }

    // here is where we write the buffer to a file
    if ((is_ready and record_all_frames and sim_is_running) or draw_this_frame) {
      static int frameno = 0;
      std::stringstream pngfn;
      pngfn << "img_" << std::setfill('0') << std::setw(5) << frameno << ".png";
      (void) saveFramePNG(pngfn.str());
      std::cout << "Wrote screenshot to " << pngfn.str() << std::endl;
      frameno++;
      draw_this_frame = false;
    }

    // draw the GUI
    ImGui::Render();

    // all done! swap buffers to the user can see
    glfwSwapBuffers(window);
  }

  // Cleanup
  std::cout << "Starting shutdown procedure" << std::endl;
  sim.reset();
  std::cout << "Quitting" << std::endl;
  ImGui_ImplGlfwGL3_Shutdown();
  glfwTerminate();

  return 0;
}

