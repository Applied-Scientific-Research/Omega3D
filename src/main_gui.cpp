/*
 * main_gui.cpp - Driver code for Omega3D + ImGui + Vc vortex particle method
 *                and boundary element method solver, GUI version
 *
 * (c)2017-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "Simulation.h"
#include "JsonHelper.h"
#include "Body.h"
#include "RenderParams.h"
#include "FeatureDraw.h"
#include "main_gui_functions.cpp"
#include "GuiHelper.h"

#include "json/json.hpp"
#include "imgui/imgui_internal.h"
#include "miniz/FrameBufferToImage.h"

#ifdef _WIN32
  // for glad
  #ifndef APIENTRY
    #define APIENTRY __stdcall
  #endif
  // for C++11 stuff that Windows can't get right
  #include <ciso646>
#endif
#include "glad.h"

//#include <GL/gl3w.h>    // This example is using gl3w to access OpenGL
// functions (because it is small). You may use glew/glad/glLoadGen/etc.
// whatever already works for you.
#include <GLFW/glfw3.h>

#include <cstdio>
#include <iostream>
#include <vector>
#include <iomanip>	// for setfill, setw


// execution starts here

int main(int argc, char const *argv[]) {

  std::cout << std::endl << "Omega3D GUI" << std::endl;
  if (VERBOSE) { std::cout << "  VERBOSE is on" << std::endl; }

  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;
  std::vector< std::unique_ptr<MeasureFeature> > mfeatures;
  FeatureDraw bdraw;
  FeatureDraw fdraw;
  FeatureDraw mdraw;
  size_t nframes = 0;
  static bool sim_is_running = false;
  static bool begin_single_step = false;

  // placeholder for command-line input file
  std::string command_line_input;

  // Set up primary OpenGL window
  glfwSetErrorCallback(error_callback);
  bool init = glfwInit();
  if (!init) {
    std::cout << "glfwInit failed" << std::endl;
    exit(-1);
  }

#if __APPLE__
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  #undef USE_OGL_COMPUTE
#elif defined (_WIN32) || defined (__linux__)
  const char* glsl_version = "#version 330 core";
#ifdef USE_OGL_COMPUTE
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
#else
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
#endif
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#else
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  #undef USE_OGL_COMPUTE
#endif
  GLFWwindow* window = glfwCreateWindow(1280, 720, "Omega3D GUI", nullptr, nullptr);
  if (!window) { 
  std::cout << "glfwCreateWindow created " << window << std::endl;
  exit(-1);
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  //gl3wInit();
  init = gladLoadGL();
  if (!init) {
    std::cout << "gladLoadGL failed" << std::endl;
    exit(-1);
  }

  // Setup ImGui binding
  ImGui::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  //glfwSetKeyCallback(keyboard_callback);

  //glfwSetWindowCloseCallback(window, window_close_callback);

#ifdef USE_OGL_COMPUTE
  // tell Simulation to generate the compute shader data structures
  sim.initGLcs();
#endif

  // Get and set some IO functions
  ImGuiIO& io = ImGui::GetIO();
  io.IniFilename = ".omega3d.ini";
  std::vector<std::string> recent_json_files;

  // Load Fonts
  // (there is a default font, this is only if you want to change it. see extra_fonts/README.txt for more details)

  //  1. (Optional) Call AddFont*** functions. If you don't call any, the default font will be loaded for you.
  //  2. Call GetTexDataAsAlpha8() or GetTexDataAsRGBA32() to build and retrieve pixels data.
  //  3. Upload the pixels data into a texture within your graphics system.
  //  4. Call SetTexID(my_tex_id); and pass the pointer/identifier to your texture. This value will be passed back to you during rendering to identify the texture.

  // increase the font size
  float fontSize = 20.0f;
  ImFontConfig config;
  //config.PixelSnapH = true;
  config.GlyphOffset = ImVec2(0.0f, 1.0f);
  config.SizePixels = fontSize;
  io.Fonts->AddFontDefault(&config);

  // try using my own font
  //io.Fonts->AddFontFromFileTTF("Roboto-Regular.ttf", 22);
  //io.Fonts->AddFontFromFileTTF("DroidSansMono.ttf", 20);
  //{
    //unsigned char* pixels;
    //int width, height;
    //io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
    //std::cout << "New font rasterized as " << width << " by " << height << std::endl;
  //}

  // a string to hold any error messages
  std::string sim_err_msg;

  // GUI and drawing parameters
  bool save_all_bdry = false;		// save Boundary Features every step
  bool save_all_flow = false;		// save Flow Features every step
  bool save_all_meas = false;		// save Measure Features every step
  bool save_all_vtus = false;		// save all collections the coming step
  bool export_vtk_this_frame = false;	// write set of vtu files with the current data
  std::vector<std::string> vtk_out_files; // list of just-output files
  bool save_all_imgs = false;		// save screenshot every step
  bool export_png_when_ready = false;	// write frame to png as soon as its done
  bool write_png_immediately = false;	// write frame to png right now
  std::string png_out_file;		// the name of the recently-written png
  bool show_stats_window = true;
  bool show_welcome_window = true;
  bool show_terminal_window = false;
  bool show_demo_window = false;
  bool show_json_input_window = false;
  bool show_file_output_window = false;
  //static bool show_origin = true;
  static bool is_viscous = false;

  // colors, modelview, and projection matrix for the render view
  RenderParams rparams;
  const bool is_ortho = false;
  std::vector<float> gl_projection;
  gl_projection.resize(16);
  std::vector<float> gl_mview;
  gl_mview.resize(16);
  if (is_ortho) compute_ortho_proj_mat(window, rparams.vsize, gl_projection);
  else compute_persp_proj_mat(window, rparams.vfov, gl_projection);
  compute_modelview_mat(rparams.vcx, rparams.vcy, rparams.vcz, rparams.rx, rparams.ry, gl_mview);

  // adjust some UI settings
  ImGuiStyle& style = ImGui::GetStyle();
  style.Colors[ImGuiCol_WindowBg]              = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TitleBg]               = ImVec4(0.27f, 0.27f, 0.54f, 1.00f);
  style.Colors[ImGuiCol_TitleBgActive]         = ImVec4(0.32f, 0.32f, 0.63f, 1.00f);

  // FE_OVERFLOW is triggered in ImGui
  // FE_INVALID is triggered in std::acos
  // FE_DIVBYZERO is triggered somewhere before Reflect.h:512, and only under x86, when switching from particle-only to BEM
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  //feenableexcept(FE_DIVBYZERO);

  // Load file names and paths of pre-stored sims
  std::vector<nlohmann::json> sims;
  std::vector<std::string> descriptions = {"Select a simulation"};
  LoadJsonSims(sims, descriptions, EXAMPLES_DIR);

  // Main loop
  std::cout << "Starting main loop" << std::endl;
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    
    //
    // Initialize simulation
    //

    if (not sim.is_initialized() and (sim_is_running || begin_single_step)) {

      std::cout << std::endl << "Initializing simulation" << std::endl;

      // initialize particle distributions
      for (auto const& ff: ffeatures) {
        if (ff->is_enabled()) {
          ElementPacket<float> newpacket = ff->init_elements(sim.get_ips());
          sim.add_elements( newpacket, active, lagrangian, ff->get_body() );
        }
      }

      // initialize solid objects
      for (auto const& bf : bfeatures) {
        if (bf->is_enabled()) {
          ElementPacket<float> newpacket = bf->init_elements(sim.get_ips());
          const move_t newMoveType = (bf->get_body() ? bodybound : fixed);
          sim.add_elements(newpacket, reactive, newMoveType, bf->get_body() );
        }
      }

      // initialize measurement features
      for (auto const& mf: mfeatures) {
        if (mf->is_enabled()) {
          ElementPacket<float> newpacket = mf->init_elements(rparams.tracer_scale*sim.get_ips());
          const move_t newMoveType = (mf->get_is_lagrangian() ? lagrangian : fixed);
          sim.add_elements(newpacket, inert, newMoveType, mf->get_body() );
        }
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
        // and reset
        sim.reset();
      }

      begin_single_step = false;
    }

    //
    // Update simulation
    //

    // get results of latest step, if it just completed
    bool is_ready = sim.test_for_new_results();

    // before we start again, write the vtu output
    if (is_ready and export_vtk_this_frame) {

      // split on which to write
      if (save_all_vtus) {
        // default is to save all collections to vtu files
        vtk_out_files = sim.write_vtk();
        // and don't do this next time
        save_all_vtus = false;

      } else if (save_all_bdry or save_all_flow or save_all_meas) {
        // only write select vtu files and don't echo
        (void) sim.write_vtk(-1, save_all_bdry, save_all_flow, save_all_meas);
      }

      // tell this routine next time around not to print
      export_vtk_this_frame = false;
    }

    // draw a notification box
    if (not vtk_out_files.empty()) {
      static int32_t vtkframect = 0;
      ++vtkframect;

      // draw the notification
      ImGui::SetNextWindowSize(ImVec2(10+fontSize*12, 10+fontSize*(2+vtk_out_files.size())));
      ImGui::SetNextWindowPos(ImVec2(0,0), 0);
      ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize;
      ImGui::Begin("Vtk written", NULL, window_flags);
      ImGui::Text("Wrote %ld file(s):", vtk_out_files.size());
      for (auto &thisfile : vtk_out_files) {
        ImGui::Text("  %s", thisfile.c_str());
      }
      ImGui::End();

      // make sure this isn't up for too long
      if (vtkframect == 90) {
        vtkframect = 0;
        vtk_out_files.clear();
      }
    }

    // see if we should start a new step
    if (is_ready and (sim_is_running || begin_single_step)) {

      if (save_all_bdry or save_all_flow or save_all_meas) export_vtk_this_frame = true;
      if (save_all_imgs) export_png_when_ready = true;

      // check flow for blow-up or dynamic errors
      sim_err_msg = sim.check_simulation();

      if (sim_err_msg.empty()) {
        // the last simulation step was fine, OK to continue
        // generate new particles from emitters
        for (auto const& ff: ffeatures) {
          if (ff->is_enabled()) {
            ElementPacket<float> newpacket = ff->step_elements(sim.get_ips());
            // echo any errors
             sim.add_elements( newpacket, active, lagrangian, ff->get_body() );
          }
        }

        for (auto const& mf: mfeatures) {
          if (mf->is_enabled()) {
            move_t newMoveType = fixed;
            if (mf->moves() or mf->emits()) {
              newMoveType = lagrangian;
            }
            sim.add_elements( mf->step_elements(rparams.tracer_scale*sim.get_ips()), inert, newMoveType, mf->get_body() );
          }
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
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiCond_FirstUseEver);
      if (ImGui::BeginPopupModal("Simulation error occurred")) {
        ImGui::Spacing();
        ImGui::TextWrapped("%s", sim_err_msg.c_str());
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
      mouse_callback(window, io, &rparams.vcx, &rparams.vcy, &rparams.vcz, &rparams.rx, &rparams.ry, &rparams.vsize);
    }

    // check for keypresses to toggle state
    //if (not io.WantCaptureKeyboard) {
      //keyboard_callback(
    //}

    //
    // The main window
    //
    {

    ImGui::SetNextWindowSize(ImVec2(140+fontSize*24,100+fontSize*12), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImVec2(20,20), ImGuiCond_FirstUseEver);
    ImGui::Begin("Omega3D");
    ImGui::Spacing();

    // Select pre-populated simulations
    {
      int currentItemIndex = 0;
      const char* currentItem = descriptions[currentItemIndex].c_str();
      static ImGuiComboFlags flags = 0;
      // The second parameter is the label previewed before opening the combo.
      if (ImGui::BeginCombo("", currentItem, flags)) {
        for (size_t n = 0; n < descriptions.size(); n++)
        {
          bool is_selected = (currentItem == descriptions[n].c_str());
          if (ImGui::Selectable(descriptions[n].c_str(), is_selected)) {
            currentItem = descriptions[n].c_str();
            currentItemIndex = n;
          }
          if (is_selected) {
              ImGui::SetItemDefaultFocus();   // Set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
          }
        }
        ImGui::EndCombo();
      }

      if (currentItemIndex) {

        // stop and clear before loading
        sim.reset();
        bfeatures.clear();
        ffeatures.clear();
        mfeatures.clear();

        // load and report
        parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, sims[currentItemIndex-1]);
        
        // clear and remake the draw geometry
        std::cout << "Loading drawing info for features..." << std::endl;
        bdraw.clear_elements();
        for (auto const& bf : bfeatures) {
          if (bf->is_enabled()) {
            bdraw.add_elements( bf->get_draw_packet(), bf->is_enabled() );
          }
        }
        fdraw.clear_elements();
        for (auto const& ff : ffeatures) {
          if (ff->is_enabled()) {
            fdraw.add_elements( ff->get_draw_packet(), ff->is_enabled() );
          }
        }
        mdraw.clear_elements();
        for (auto const& mf : mfeatures) {
          if (mf->is_enabled()) {
            mdraw.add_elements( mf->get_draw_packet(), mf->is_enabled() );
          }
        }
        // finish setting up and run
        is_viscous = sim.get_diffuse();
        currentItemIndex = 0;
      }
    }

    // or load a simulation from a JSON file
    ImGui::SameLine();
    if (ImGui::Button("Load from json", ImVec2(10+fontSize*8,0))) show_json_input_window = true;

    if (show_json_input_window) {

      // FileIO is now a modal, see code in lib/imgui/ImguiWindowsFileIO.cpp
      ImGui::OpenPopup("FileIO");

      bool try_it = false;
      static std::string infile = "input.json";
      if (fileIOWindow(try_it, infile, recent_json_files, "Open", {"*.json", "*.*"}, true, ImVec2(200+26*fontSize,300))) {
        show_json_input_window = false;

        if (try_it and !infile.empty()) {
          // remember
          recent_json_files.push_back( infile );

          // stop and clear before loading
          sim.reset();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();

          // load and report
          nlohmann::json j = read_json(infile);
          parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, j);

          // clear and remake the draw geometry
          std::cout << "Loading drawing info for features..." << std::endl;
          bdraw.clear_elements();
          for (auto const& bf : bfeatures) {
            if (bf->is_enabled()) {
              bdraw.add_elements( bf->get_draw_packet(), bf->is_enabled() );
            }
          }
          fdraw.clear_elements();
          for (auto const& ff : ffeatures) {
            if (ff->is_enabled()) {
              fdraw.add_elements( ff->get_draw_packet(), ff->is_enabled() );
            }
          }
          mdraw.clear_elements();
          for (auto const& mf : mfeatures) {
            if (mf->is_enabled()) {
              mdraw.add_elements( mf->get_draw_packet(), mf->is_enabled() );
            }
          }
        
          // we have to manually set this variable
          is_viscous = sim.get_diffuse();

          // autostart if file requests it
          if (sim.autostart()) sim_is_running = true;

          // check and possibly resize the window to match the saved resolution
          resize_to_resolution(window, rparams.width, rparams.height);
        }
      }
    }
    ImGui::Spacing();

    // or load the sim from the command-line (do this once)
    if (argc == 2 and command_line_input.empty()) {

      // stop and clear before loading
      sim.reset();
      bfeatures.clear();
      ffeatures.clear();
      mfeatures.clear();

      // load and report
      command_line_input = argv[1];
      nlohmann::json j = read_json(command_line_input);
      parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, j);

      // clear and remake the draw geometry
      std::cout << "Loading drawing info for features..." << std::endl;
      bdraw.clear_elements();
      for (auto const& bf : bfeatures) {
        if (bf->is_enabled()) {
          bdraw.add_elements( bf->get_draw_packet(), bf->is_enabled() );
        }
      }
      fdraw.clear_elements();
      for (auto const& ff : ffeatures) {
        if (ff->is_enabled()) {
          fdraw.add_elements( ff->get_draw_packet(), ff->is_enabled() );
        }
      }
      mdraw.clear_elements();
      for (auto const& mf : mfeatures) {
        if (mf->is_enabled()) {
          mdraw.add_elements( mf->get_draw_packet(), mf->is_enabled() );
        }
      }

      // we have to manually set this variable
      is_viscous = sim.get_diffuse();

      // autostart if file requests it
      if (sim.autostart()) sim_is_running = true;

      // check and possibly resize the window to match the saved resolution
      resize_to_resolution(window, rparams.width, rparams.height);

      // we don't need the welcome banner
      show_welcome_window = false;
    }


    //if (ImGui::CollapsingHeader("Simulation globals", ImGuiTreeNodeFlags_DefaultOpen))
    if (ImGui::CollapsingHeader("Simulation globals")) {

      // save current versions, so we know which changed
      const float current_re = sim.get_re();
      const float current_dt = sim.get_dt();

      ImGui::Checkbox("Fluid is viscous (diffuses)", &is_viscous);
      ImGui::SameLine();
      ShowHelpMarker("If checked, simulation will add particles to solve diffusion equation; if unchecked, simulation will run fast, but quickly lose accuracy.");

      if (sim.is_initialized() and is_viscous) {
        ImGui::Text("** Until a reset, Time step and Reynolds number move together **");
      }

      // set input widget width for this and the next few
      ImGui::PushItemWidth(-240);

      ImGui::SliderFloat("Time step", sim.addr_dt(), 0.001f, 0.1f, "%.4f", 2.0f);
      ImGui::SameLine();
      ShowHelpMarker("Adjust how far into the future each step must simulate. Smaller means better accuracy but slower, larger means lower accuracy but faster.");

      if (is_viscous) {
        sim.set_diffuse(true);

        // and let user choose Reynolds number
        ImGui::SliderFloat("Reynolds number", sim.addr_re(), 10.0f, 2000.0f, "%.1f", 2.0f);
        ImGui::SameLine();
        ShowHelpMarker("Reynolds number is the inverse of viscosity; so larger means less viscosity, smaller particles, and longer run time, but more detail.");

        // if Reynolds number or time step change during a run, adjust the other to keep particle spacing constant
        if (sim.is_initialized()) {
          if (current_re != sim.get_re()) {
            // change dt
            *(sim.addr_dt()) = current_dt * sim.get_re() / current_re;
          } else if (current_dt != sim.get_dt()) {
            // change Re
            *(sim.addr_re()) = current_re * sim.get_dt() / current_dt;
          }
        }

        ImGui::Text("Particle spacing for these settings is %g", sim.get_ips());

      } else {
        static float my_ips = 0.03141;
        ImGui::SliderFloat("Particle spacing", &my_ips, 0.001f, 0.1f, "%.3f", 2.0f);
        ImGui::SameLine();
        ShowHelpMarker("Sets the average size and distance between particles.");
        // change underlying Re when this changes
        sim.set_re_for_ips(my_ips);
        my_ips = sim.get_ips();
      }
      ImGui::InputFloat3("Flow vector", sim.addr_fs());
      ImGui::SameLine();
      ShowHelpMarker("Also called freestream, this is a uniform wind blowing along this vector.");

      ImGui::PopItemWidth();

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
    } // end Simulation Globals

    ImGui::Spacing();
    //if (ImGui::CollapsingHeader("Flow structures", ImGuiTreeNodeFlags_DefaultOpen)) 
    if (ImGui::CollapsingHeader("Startup structures")) {

      if (ffeatures.size() + bfeatures.size() == 0) {
        ImGui::Text("Add flow or boundry features (like vortex blobs and solid objects) here, then click RUN.");
      }

      ImGui::Spacing();

      // button and window for adding new boundary objects
      static bool create_bdry_f = false;
      if (ImGui::Button("Add boundary")) { create_bdry_f = true; }
      if (create_bdry_f) {
        ImGui::SetNextWindowPos(ImVec2(100, 100), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(400, 275), ImGuiCond_FirstUseEver);
        ImGui::Begin("New boundary structure");
        int status = BoundaryFeature::draw_creation_gui(bfeatures, sim);
        if (status == 1) {
            bdraw.add_elements(bfeatures.back()->get_draw_packet(), bfeatures.back()->is_enabled());
            create_bdry_f = false;
        } else if (status == 2) {
            create_bdry_f = false;
        }
        ImGui::End();
      }

      // button and window for adding new flow structures
      ImGui::SameLine();
      static bool create_flow_f = false;
      if (ImGui::Button("Add vortex")) { create_flow_f = true; }
      if (create_flow_f) {
        ImGui::SetNextWindowPos(ImVec2(100, 100), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(500, 300), ImGuiCond_FirstUseEver);
        ImGui::Begin("New flow structure");
        int status = FlowFeature::draw_creation_gui(ffeatures, sim.get_ips());
        if (status == 1) {
            fdraw.add_elements(ffeatures.back()->get_draw_packet(), ffeatures.back()->is_enabled());
            create_flow_f = false;
        } else if (status == 2) {
            create_flow_f = false;
        }
        ImGui::End();
      }

      // button and window for adding new measurement objects
      ImGui::SameLine();
      static bool create_ms_f = false;
      if (ImGui::Button("Add measurement")) { create_ms_f = true; }
      if (create_ms_f) {
        ImGui::SetNextWindowPos(ImVec2(100, 100), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(400, 200), ImGuiCond_FirstUseEver);
        ImGui::Begin("New measurement structure");
        int status = MeasureFeature::draw_creation_gui(mfeatures, sim.get_ips(), rparams.tracer_scale);
        if (status == 1) {
          mdraw.add_elements(mfeatures.back()->get_draw_packet(), mfeatures.back()->is_enabled());
          create_ms_f = false;
        } else if (status == 2) {
          create_ms_f = false;
        }
        ImGui::End();
      }

      ImGui::Spacing();

      // list existing flow features here
      static std::unique_ptr<FlowFeature> tmpff = nullptr;
      static int edit_feat_index = -1;
      int del_feat_index = -1;
      bool redraw = false;
      int buttonIDs = 10;
      FlowFeature::draw_feature_list(ffeatures, tmpff, edit_feat_index, del_feat_index, redraw, buttonIDs);
      
      if (tmpff) {
        ImGui::OpenPopup("Edit flow feature");
        ImGui::SetNextWindowSize(ImVec2(400,275), ImGuiCond_FirstUseEver);
        if (ImGui::BeginPopupModal("Edit flow feature")) {
          bool fin = false;
          if (tmpff->draw_info_gui("Edit", sim.get_ips())) {
            tmpff->generate_draw_geom();
            ffeatures[edit_feat_index].reset(nullptr);
            ffeatures[edit_feat_index] = std::move(tmpff);
            redraw = true;
            fin = true;
          }
          ImGui::SameLine();
          if (ImGui::Button("Cancel", ImVec2(120,0))) {
            fin = true;
          }
          if (fin) {
            edit_feat_index = -1;
            tmpff = nullptr;
            ImGui::CloseCurrentPopup();
          }
        ImGui::EndPopup();
        }
      }

      if (del_feat_index > -1) {
        std::cout << "Asked to delete flow feature " << del_feat_index << std::endl;
        ffeatures.erase(ffeatures.begin()+del_feat_index);
        redraw = true;
        del_feat_index = -1;
      }

      if (redraw) {
        fdraw.clear_elements();
        for (auto const& ff : ffeatures) {
          if (ff->is_enabled()) {
            fdraw.add_elements( ff->get_draw_packet(), ff->is_enabled() );
          }
        }
        redraw = false;
      }

      // list existing boundary features here
      static std::unique_ptr<BoundaryFeature> tmpbf = nullptr;
      BoundaryFeature::draw_feature_list(bfeatures, tmpbf, edit_feat_index, del_feat_index, redraw, buttonIDs);

      if (tmpbf) {
        ImGui::OpenPopup("Edit boundary feature");
        ImGui::SetNextWindowSize(ImVec2(400,275), ImGuiCond_FirstUseEver);
        if (ImGui::BeginPopupModal("Edit boundary feature")) {
          bool fin = false;
          // Currently cannot edit body. This will require rethinking on how we manage the Body Class.
          if (tmpbf->draw_info_gui("Edit")) {
            tmpbf->create();
            tmpbf->generate_draw_geom();
            bfeatures[edit_feat_index].reset();
            bfeatures[edit_feat_index] = std::move(tmpbf);
            redraw = true;
            fin = true;
          }
          ImGui::SameLine();
          if (ImGui::Button("Cancel", ImVec2(120,0))) { fin = true; }
          if (fin) {
            edit_feat_index = -1;
            tmpbf = nullptr;
            ImGui::CloseCurrentPopup();
          }
          ImGui::EndPopup();
        }
      }

      if (del_feat_index > -1) {
        std::cout << "Asked to delete boundary feature " << del_feat_index<< std::endl;
        bfeatures.erase(bfeatures.begin()+del_feat_index);
        del_feat_index = -1;
        redraw = true;
      }

      if (redraw) {
        // clear out and re-make all boundary draw geometry
        bdraw.clear_elements();
        for (auto const& bf : bfeatures) {
          if (bf->is_enabled()) { 
            bdraw.add_elements( bf->get_draw_packet(), bf->is_enabled() );
          }
        }
        redraw = false;
      }

      // list existing measurement features here
      static std::unique_ptr<MeasureFeature> tmpmf = nullptr;
      MeasureFeature::draw_feature_list(mfeatures, tmpmf, edit_feat_index, del_feat_index, redraw, buttonIDs);
   
      if (tmpmf) {
        ImGui::OpenPopup("Edit measure feature");
        ImGui::SetNextWindowSize(ImVec2(400,275), ImGuiCond_FirstUseEver);
        if (ImGui::BeginPopupModal("Edit measure feature")) {
          bool fin = false;
          if (tmpmf->draw_info_gui("Edit", rparams.tracer_scale, sim.get_ips())) {
            tmpmf->generate_draw_geom();
            mfeatures[edit_feat_index].reset();
            mfeatures[edit_feat_index] = std::move(tmpmf);
            redraw = true;
            fin = true;
          }
          ImGui::SameLine();
          if (ImGui::Button("Cancel", ImVec2(120,0))) { fin = true; }
          if (fin) {
            edit_feat_index = -1;
            tmpmf = nullptr;
            ImGui::CloseCurrentPopup();
          }
        ImGui::EndPopup();
        }
      }

      if (del_feat_index > -1) {
        std::cout << "Asked to delete measurement feature " << del_feat_index<< std::endl;
        mfeatures.erase(mfeatures.begin()+del_feat_index);
        redraw = true;
        del_feat_index = -1;
      }

      if (redraw) {
        mdraw.clear_elements();
        for (auto const& mf : mfeatures) {
          if (mf->is_enabled()) {
            mdraw.add_elements( mf->get_draw_packet(), mf->is_enabled() );
          }
        }
        redraw = false;
      }
    } // end structure entry


    // Rendering parameters, under a header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Rendering controls")) { draw_render_gui(rparams); }

    // Solver parameters, under its own header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Solver parameters (advanced)")) { sim.draw_advanced(); }

    // Output buttons, under a header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Save output")) {

      // save setup
      ImGui::Spacing();
      ImGui::Text("Save simulation setup:");
      ImGui::SameLine();
      if (ImGui::Button("to JSON", ImVec2(10+4*fontSize,0))) show_file_output_window = true;

      // save current data
      ImGui::Separator();
      ImGui::Spacing();
      ImGui::Text("Save current step:");

      ImGui::SameLine();
      if (ImGui::Button("All to VTU", ImVec2(10+7*fontSize,0))) {
        save_all_vtus = true;
        export_vtk_this_frame = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Screenshot to PNG", ImVec2(10+10*fontSize,0))) write_png_immediately = true;

      // save data regularly
      ImGui::Separator();
      ImGui::Spacing();
      ImGui::Text("Save every step:");

      ImGui::Indent();
      ImGui::Checkbox("Boundary features (to VTU)", &save_all_bdry);
      ImGui::Checkbox("Flow features (to VTU)", &save_all_flow);
      ImGui::Checkbox("Measure features (to VTU)", &save_all_meas);
      ImGui::Checkbox("Window screenshot (to PNG)", &save_all_imgs);
      ImGui::Unindent();
    }

    if (show_file_output_window) {
      bool try_it = false;
      static std::string outfile = "file_name.json";

      // FileIO is now a modal, see code in lib/imgui/ImguiWindowsFileIO.cpp
      ImGui::OpenPopup("FileIO");

      if (fileIOWindow(try_it, outfile, recent_json_files, "Save", {"*.json"}, false, ImVec2(200+26*fontSize,300))) {
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
    nframes++;

    // check vs. end conditions, if present
    if (sim.test_vs_stop_async()) {
      // just pause sim
      sim_is_running = false;

      // fully quit if asked
      if (sim.quitonstop()) {
        if (is_ready) {
          // this simulation step has finished, write png and exit
          write_png_immediately = true;

          // tell glfw to close the window next time around
          glfwSetWindowShouldClose(window, GLFW_TRUE);
        }
      }
    }

    // all the other stuff
    {
      ImGui::Spacing();
      ImGui::Separator();
      ImGui::Spacing();

      if (sim_is_running) {
        //ImGui::Text("Simulation is running...step = %ld, time = %g", sim.get_nstep(), sim.get_time());
        if (ImGui::Button("PAUSE", ImVec2(200,20+fontSize))) sim_is_running = false;
        // space bar pauses
        if (ImGui::IsKeyPressed(32) and not show_file_output_window) sim_is_running = false;
      } else if (!(sim_is_running) && !(sim.test_for_new_results())) {
        ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
        ImGui::Button("FINISHING STEP", ImVec2(200, 20+fontSize));
        ImGui::PopItemFlag();
      } else {
        //ImGui::Text("Simulation is not running, step = %ld, time = %g", sim.get_nstep(), sim.get_time());
        if (ImGui::Button("RUN", ImVec2(200,20+fontSize))) sim_is_running = true;
        ImGui::SameLine();
        if (ImGui::Button("Step", ImVec2(120,0))) begin_single_step = true;
        // and space bar resumes
        if (ImGui::IsKeyPressed(32) and not show_file_output_window) sim_is_running = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Reset", ImVec2(120,0))) {
        std::cout << std::endl << "Reset requested" << std::endl;
        // remove all particles and reset timer
        sim.reset();
        std::cout << "Reset complete" << std::endl;
      }

      //ImGui::Spacing();
      //ImGui::Separator();
      /*
      ImGui::Text("Open additional windows");
      if (ImGui::Button("Plot statistics")) show_stats_window ^= 1;
      ImGui::SameLine();
      if (ImGui::Button("Show terminal output")) show_terminal_window ^= 1;
      ImGui::SameLine();
      */

      //if (ImGui::Button("ImGui Samples")) show_demo_window ^= 1;
      // use ASCII table for number: http://www.asciitable.com/
      // but use CAPITAL letter for a letter, jesus, really?!?
      //if (ImGui::IsKeyPressed(84) and not show_file_output_window) show_demo_window ^= 1;

      //ImGui::Text("Draw frame rate: %.2f ms/frame (%.1f FPS)",
      //            1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

      //ImGui::Text("Number of panels: %ld  particles: %ld  field points: %ld",
      //            sim.get_npanels(), sim.get_nparts(), sim.get_nfldpts());
    }


    // done drawing the UI window
    ImGui::End();
    }

    // Rendering
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(rparams.clear_color[0], rparams.clear_color[1], rparams.clear_color[2], rparams.clear_color[3]);
    glClear(GL_COLOR_BUFFER_BIT);

    // Show the welcome window
    if (show_welcome_window) { show_welcome_window = draw_welcome_window(io.DisplaySize.x, io.DisplaySize.y); }

    // Show the simulation stats in the corner
    //if (nframes < 10){ std::cout << "show_stats_window: " << show_stats_window << std::endl; }
    if (show_stats_window) { draw_stats_window(sim.get_npanels(), sim.get_nfldpts(), sim.get_nstep(), 
                                               sim.get_time(), sim.get_nparts(), &show_stats_window,
                                               fontSize, display_h); }

    // Show the terminal output of the program
    if (show_terminal_window) {
      ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiCond_FirstUseEver);
      ImGui::Begin("Terminal", &show_terminal_window);
      ImGui::Text("Hello");
      ImGui::End();
    }

    // Show the ImGui test window. Most of the sample code is in ImGui::ShowTestWindow()
    if (show_demo_window) {
      ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
      ImGui::ShowDemoWindow();
    }

#ifdef USE_OGL_COMPUTE
    // use compute shaders to advance the simulation
    sim.computeGL();
#endif

    // draw the simulation: panels and particles
    if (is_ortho) compute_ortho_proj_mat(window, rparams.vsize, gl_projection);
    else compute_persp_proj_mat(window, rparams.vfov, gl_projection);
    compute_modelview_mat(rparams.vcx, rparams.vcy, rparams.vcz, rparams.rx, rparams.ry, gl_mview);
    sim.drawGL(gl_mview, gl_projection, rparams);

    // if simulation has not been initted, draw the features instead!
    if (not sim.is_initialized()) {
      // append draw geometries to FeatureDraw object
      for (auto const& bf : bfeatures) {
        // Whatever happens here should happen with measure/flow features
        if (bf->is_enabled()) {
          // what should we do differently?
        }
      }

      // and draw
      //std::cout << "Main_Gui pre-sim draw call" << std::endl;
      bdraw.drawGL(gl_mview, gl_projection, rparams, true);
      fdraw.drawGL(gl_mview, gl_projection, rparams, false);
      mdraw.drawGL(gl_mview, gl_projection, rparams, true);
    }

    // here is where we write the buffer to a file
    if ((is_ready and export_png_when_ready) or write_png_immediately) {
      static int frameno = 0;
      std::stringstream pngfn;
      pngfn << "img_" << std::setfill('0') << std::setw(5) << frameno << ".png";
      png_out_file = pngfn.str();
      (void) saveFramePNG(png_out_file);
      std::cout << "Wrote screenshot to " << png_out_file << std::endl;
      frameno++;
      // no need to tell the user every frame
      if (export_png_when_ready) png_out_file.clear();
      write_png_immediately = false;
      export_png_when_ready = false;
    }

    // if we're just drawing this one frame, then announce that we wrote it
    if (not png_out_file.empty()) {
      static int32_t pngframect = 0;
      ++pngframect;

      // draw the notification
      ImGui::SetNextWindowSize(ImVec2(10+fontSize*12, 10+fontSize*2));
      ImGui::SetNextWindowPos(ImVec2(0,0), 0);
      ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize;
      ImGui::Begin("Png written", NULL, window_flags);
      ImGui::Text("Wrote %s", png_out_file.c_str());
      ImGui::End();

      // make sure this isn't up for too long
      if (pngframect == 90) {
        pngframect = 0;
        png_out_file.clear();
      }
    }

    // draw the GUI
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    // all done! swap buffers to the user can see
    glfwMakeContextCurrent(window);
    glfwSwapBuffers(window);
  }

  // Cleanup
  std::cout << "Starting shutdown procedure" << std::endl;
  sim.reset();
  std::cout << "Quitting" << std::endl;
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}

