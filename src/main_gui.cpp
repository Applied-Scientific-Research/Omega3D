/*
 * main_gui.cpp - Driver code for Omega3D + ImGui + Vc vortex particle method
 *                and boundary element method solver, GUI version
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "Simulation.h"
#include "JsonHelper.h"
#include "Body.h"
#include "RenderParams.h"
#include "json/json.hpp"
#include "main_gui_functions.cpp"

#ifdef _WIN32
  // for glad
  #ifndef APIENTRY
    #define APIENTRY __stdcall
  #endif
  // for C++11 stuff that Windows can't get right
  #include <ciso646>
#endif
#include "glad.h"

// header-only immediate-mode GUI
#include "GuiHelper.h"

// header-only png writing
#include "stb/FrameBufferToImage.h"

//#include <GL/gl3w.h>    // This example is using gl3w to access OpenGL
// functions (because it is small). You may use glew/glad/glLoadGen/etc.
// whatever already works for you.
#include <GLFW/glfw3.h>

#include <cstdio>
#include <iostream>
#include <vector>
#include <iomanip>	// for setfill, setw

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
#elif _WIN32
  const char* glsl_version = "#version 330 core";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#else
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
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
  //sim.set_status_file_name("status.dat");

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
  bool export_vtk_this_frame = false;	// write a vtk with the current data
  std::vector<std::string> vtk_out_files; // list of just-output files
  bool draw_this_frame = false;		// draw the frame as soon as its done
  std::string png_out_file;		// the name of the recently-written png
  bool record_all_frames = false;	// save a frame when a new one is ready
  bool show_stats_window = true;
  bool show_welcome_window = true;
  bool show_terminal_window = false;
  bool show_demo_window = false;
  bool show_json_input_window = false;
  bool show_file_output_window = false;
  //static bool show_origin = true;
  static bool is_viscous = false;

  // colors and projection matrix for the render view
  RenderParams rparams;
  std::vector<float> gl_projection;
  compute_ortho_proj_mat(window, rparams.vcx, rparams.vcy, &rparams.vsize, gl_projection);

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
        if (ff->is_enabled()) sim.add_particles( ff->init_particles(sim.get_ips()) );
      }

      // initialize solid objects
      for (auto const& bf : bfeatures) {
        if (bf->is_enabled()) sim.add_boundary( bf->get_body(), bf->init_elements(sim.get_ips()) );
      }

      // initialize measurement features
      for (auto const& mf: mfeatures) {
        if (mf->is_enabled()) sim.add_fldpts( mf->init_particles(rparams.tracer_scale*sim.get_ips()), mf->moves() );
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
      vtk_out_files = sim.write_vtk();
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

      // check flow for blow-up or dynamic errors
      sim_err_msg = sim.check_simulation();

      if (sim_err_msg.empty()) {
        // the last simulation step was fine, OK to continue

        // generate new particles from emitters
        for (auto const& ff : ffeatures) {
          if (ff->is_enabled()) sim.add_particles( ff->step_particles(sim.get_ips()) );
        }
        for (auto const& mf : mfeatures) {
          if (mf->is_enabled()) sim.add_fldpts( mf->step_particles(rparams.tracer_scale*sim.get_ips()), true );
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
      if (ImGui::BeginCombo("", currentItem, flags)) // The second parameter is the label previewed before opening the combo.
      {
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

      if(currentItemIndex) {
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, sims[currentItemIndex-1]);
          is_viscous = sim.get_diffuse();
          currentItemIndex = 0;
          sim_is_running = true;
      }
    }

    // or load a simulation from a JSON file
    ImGui::SameLine();
    if (ImGui::Button("Load from json", ImVec2(10+fontSize*8,0))) show_json_input_window = true;

    if (show_json_input_window) {
      bool try_it = false;
      static std::string infile = "input.json";

      if (fileIOWindow( try_it, infile, recent_json_files, "Open", {"*.json", "*.*"}, true, ImVec2(200+26*fontSize,300))) {
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

          // we have to manually set this variable
          is_viscous = sim.get_diffuse();

          // run one step so we know what we have, or autostart
          if (sim.autostart()) {
            sim_is_running = true;
          } else {
            begin_single_step = true;
          }

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

      command_line_input = argv[1];
      nlohmann::json j = read_json(command_line_input);
      parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, j);

      // we have to manually set this variable
      is_viscous = sim.get_diffuse();

      // run one step so we know what we have, or autostart
      if (sim.autostart()) {
        sim_is_running = true;
      } else {
        begin_single_step = true;
      }

      // check and possibly resize the window to match the saved resolution
      resize_to_resolution(window, rparams.width, rparams.height);

      // we don't need the welcome banner
      show_welcome_window = false;
    }


    //if (ImGui::CollapsingHeader("Simulation globals", ImGuiTreeNodeFlags_DefaultOpen)) {
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
    //if (ImGui::CollapsingHeader("Flow structures", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (ImGui::CollapsingHeader("Startup structures")) {

      if (ffeatures.size() + bfeatures.size() == 0) {
        ImGui::Text("Add flow or boundry features (like vortex blobs and solid objects) here, then click RUN.");
      }

      ImGui::Spacing();

      // button and modal window for adding new boundary objects
      if (ImGui::Button("Add boundary")) ImGui::OpenPopup("New boundary structure");
      ImGui::SetNextWindowSize(ImVec2(400,275), ImGuiCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New boundary structure")) {
        if (BoundaryFeature::draw_creation_gui(bfeatures, sim)) {
          // Eventually draw figure pre-sim
        }
        ImGui::EndPopup();
      }

      // button and modal window for adding new flow structures
      ImGui::SameLine();
      if (ImGui::Button("Add vortex")) ImGui::OpenPopup("New flow structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiCond_FirstUseEver);
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
          case 0: {
            // a blob of multiple vortons
            ImGui::InputFloat3("strength", vstr);
            ImGui::SliderFloat("radius", &rad, sim.get_ips(), 10.0f*sim.get_ips(), "%.4f");
            ImGui::SliderFloat("softness", &soft, sim.get_ips(), rad, "%.4f");
            guess_n = 4.1888f * std::pow((2.0f*rad+soft)/sim.get_ips(), 3);
            ImGui::Spacing();
            ImGui::TextWrapped("This feature will add about %d particles", guess_n);
            ImGui::Spacing();
            if (ImGui::Button("Add vortex blob")) {
              ffeatures.emplace_back(std::make_unique<VortexBlob>(xc[0],xc[1],xc[2], vstr[0],vstr[1],vstr[2], rad, soft));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            // it would be nice to be able to put this all in
            //SingleParticle::draw_creation_gui();
            } break;

          case 1: {
            // random particles in a block
            ImGui::SliderInt("number", &npart, 10, 100000);
            ImGui::SliderFloat3("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
            ImGui::SliderFloat("strength magnitude", &strmag, 0.01f, 10.0f, "%.3f", 2.0f);
            ImGui::Spacing();
            ImGui::TextWrapped("This feature will add %d particles", npart);
            ImGui::Spacing();
            if (ImGui::Button("Add random vorticies")) {
              ffeatures.emplace_back(std::make_unique<BlockOfRandom>(xc[0],xc[1],xc[2], xs[0],xs[1],xs[2], strmag, npart));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 2: {
            // oriented singular vortex ring
            ImGui::InputFloat3("direction", vstr);
            ImGui::SliderFloat("circulation", &circ, 0.001f, 10.0f, "%.3f");
            ImGui::SliderFloat("radius", &rad, 3.0f*sim.get_ips(), 10.0f, "%.3f");
            guess_n = 1 + (2.0f * 3.1416f * rad / sim.get_ips());
            ImGui::Spacing();
            ImGui::TextWrapped("This feature will add about %d particles", guess_n);
            ImGui::Spacing();
            if (ImGui::Button("Add singular vortex ring")) {
              ffeatures.emplace_back(std::make_unique<SingularRing>(xc[0],xc[1],xc[2], vstr[0],vstr[1],vstr[2], rad, circ));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 3: {
            // oriented thick-cored vortex ring
            ImGui::InputFloat3("direction", vstr);
            ImGui::SliderFloat("circulation", &circ, 0.001f, 10.0f, "%.4f");
            ImGui::SliderFloat("radius", &rad, 3.0f*sim.get_ips(), 10.0f, "%.3f");
            ImGui::SliderFloat("thickness", &soft, sim.get_ips(), 10.0f*sim.get_ips(), "%.4f");
            guess_n = (1 + (2.0f * 3.1416f * rad / sim.get_ips())) * std::pow(soft / sim.get_ips(), 2);
            ImGui::Spacing();
            ImGui::TextWrapped("This feature will add about %d particles", guess_n);
            ImGui::Spacing();
            if (ImGui::Button("Add thick vortex ring")) {
              ffeatures.emplace_back(std::make_unique<ThickRing>(xc[0],xc[1],xc[2], vstr[0],vstr[1],vstr[2], rad, soft, circ));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      } // end popup new flow structures


      // button and modal window for adding new measurement objects
      ImGui::SameLine();
      if (ImGui::Button("Add measurement")) ImGui::OpenPopup("New measurement structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New measurement structure"))
      {
        static int item = 0;
        const char* items[] = { "single point/tracer", "streakline", "circle of tracers", "line of tracers", "measurement line", "2D grid" };
        ImGui::Combo("type", &item, items, 6);

        static float xc[3] = {0.0f, 0.0f, 0.0f};
        static float xf[3] = {0.0f, 1.0f, 0.0f};
        static bool is_lagrangian = true;
        static float rad = 2.0 * sim.get_ips();

        // show different inputs based on what is selected
        switch(item) {
          case 0: {
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
            } break;
          case 1: {
            // a tracer emitter
            ImGui::InputFloat3("position", xc);
            ImGui::TextWrapped("This feature will add 1 tracer emitter");
            if (ImGui::Button("Add streakline")) {
              mfeatures.emplace_back(std::make_unique<TracerEmitter>(xc[0], xc[1], xc[2]));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
          case 2: {
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
            } break;
          case 3: {
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
            } break;
          case 4: {
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
            } break;
          case 5: {
            // a static, 2d measurement grid
            static float xg[3] = {1.0f, 0.0f, 0.0f};
            static float dx[2] = {0.1f, 0.1f};
            ImGui::InputFloat3("corner", xc);
            ImGui::InputFloat3("axis 1", xf);
            ImGui::InputFloat3("axis 2", xg);
            ImGui::InputFloat2("dx", dx);
            ImGui::TextWrapped("This feature will add about %d field points",
                               1+(int)(std::sqrt(xf[0]*xf[0]+xf[1]*xf[1]+xf[2]*xf[2])*
                                       std::sqrt(xg[0]*xg[0]+xg[1]*xg[1]+xg[2]*xg[2])/
                                       (dx[0]*dx[1])));
            if (ImGui::Button("Add 2D grid of measurement points")) {
              mfeatures.emplace_back(std::make_unique<Grid2dPoints>(xc[0], xc[1], xc[2], xf[0], xf[1], xf[2], xg[0], xg[1], xg[2], dx[0], dx[1]));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
        } // end measurement structures 

      ImGui::Spacing();
      int buttonIDs = 10;

      // list existing flow features here
      int del_this_item = -1;
      for (int i=0; i<(int)ffeatures.size(); ++i) {

        ImGui::PushID(++buttonIDs);
        ImGui::Checkbox("", ffeatures[i]->addr_enabled());
        ImGui::PopID();
        if (ffeatures[i]->is_enabled()) {
          ImGui::SameLine();
          ImGui::Text("%s", ffeatures[i]->to_string().c_str());
        } else {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1.0f), "%s", ffeatures[i]->to_string().c_str());
        }

        // add a "remove" button at the end of the line (so it's not easy to accidentally hit)
        ImGui::SameLine();
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_item = i;
        ImGui::PopID();

        //ImGui::SameLine();
        //ImGui::PushID(++buttonIDs);
        //if (ImGui::SmallButton("edit", ImVec2(60,0))) edit_this_item = i;
        //ImGui::PopID();
      }
      if (del_this_item > -1) {
        std::cout << "Asked to delete flow feature " << del_this_item << std::endl;
        ffeatures.erase(ffeatures.begin()+del_this_item);
      }

      // list existing boundary features here
      int del_this_bdry = -1;
      for (int i=0; i<(int)bfeatures.size(); ++i) {

        ImGui::PushID(++buttonIDs);
        ImGui::Checkbox("", bfeatures[i]->addr_enabled());
        ImGui::PopID();
        if (bfeatures[i]->is_enabled()) {
          ImGui::SameLine();
          ImGui::Text("%s", bfeatures[i]->to_string().c_str());
        } else {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1.0f), "%s", bfeatures[i]->to_string().c_str());
        }

        // add a "remove" button at the end of the line (so it's not easy to accidentally hit)
        ImGui::SameLine();
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_bdry = i;
        ImGui::PopID();
      }
      if (del_this_bdry > -1) {
        std::cout << "Asked to delete boundary feature " << del_this_bdry << std::endl;
        bfeatures.erase(bfeatures.begin()+del_this_bdry);
      }

      // list existing measurement features here
      int del_this_measure = -1;
      for (int i=0; i<(int)mfeatures.size(); ++i) {

        ImGui::PushID(++buttonIDs);
        ImGui::Checkbox("", mfeatures[i]->addr_enabled());
        ImGui::PopID();
        if (mfeatures[i]->is_enabled()) {
          ImGui::SameLine();
          ImGui::Text("%s", mfeatures[i]->to_string().c_str());
        } else {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1.0f), "%s", mfeatures[i]->to_string().c_str());
        }

        // add a "remove" button at the end of the line (so it's not easy to accidentally hit)
        ImGui::SameLine();
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_measure = i;
        ImGui::PopID();
      }
      if (del_this_measure > -1) {
        std::cout << "Asked to delete measurement feature " << del_this_measure << std::endl;
        mfeatures.erase(mfeatures.begin()+del_this_measure);
      }

      //if (ffeatures.size() + bfeatures.size() + mfeatures.size() == 0) {
      //  ImGui::Text("none");
      //}

    } // end structure entry


    // Rendering parameters, under a header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Rendering controls")) {
      ImGui::ColorEdit3("positive circulation", rparams.pos_circ_color);
      ImGui::ColorEdit3("negative circulation", rparams.neg_circ_color);
      ImGui::ColorEdit3("feature color",        rparams.default_color);
      ImGui::ColorEdit3("background color",     rparams.clear_color);
      //ImGui::Checkbox("show origin", &show_origin);
      ImGui::SliderFloat("particle brightness", &(rparams.circ_density), 0.0001f, 1.0f, "%.4f", 4.0f);
      ImGui::SliderFloat("particle scale", &(rparams.vorton_scale), 0.01f, 1.5f, "%.2f", 2.0f);

      if (ImGui::Button("Recenter")) {
        // put everything back to center
        rparams.vcx = -0.5f;
        rparams.vcy = 0.0f;
        rparams.vsize = 2.0f;
      }

      // add button to recenter on all vorticity?
    }

    // Solver parameters, under its own header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Solver parameters (advanced)")) { sim.draw_advanced(); }

    // Output buttons, under a header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Save output")) {

      // save the simulation to a JSON or VTK file
      ImGui::Spacing();
      if (ImGui::Button("Save setup to JSON", ImVec2(20+12*fontSize,0))) show_file_output_window = true;
      ImGui::SameLine();
      // PNG output of the render frame
      if (ImGui::Button("Save screenshot to PNG", ImVec2(20+12*fontSize,0))) draw_this_frame = true;

      // next line: VTK output and record
      if (ImGui::Button("Save parts to VTU", ImVec2(20+12*fontSize,0))) export_vtk_this_frame = true;
      ImGui::SameLine();
      if (record_all_frames) {
        if (ImGui::Button("STOP RECORDING", ImVec2(20+12*fontSize,0))) {
          record_all_frames = false;
          sim_is_running = false;
        }
      } else {
        if (ImGui::Button("RECORD to PNG", ImVec2(20+12*fontSize,0))) {
          record_all_frames = true;
          sim_is_running = true;
        }
      }
    }

    if (show_file_output_window) {
      bool try_it = false;
      static std::string outfile = "output.json";

      if (fileIOWindow( try_it, outfile, recent_json_files, "Save", {"*.json", "*.*"}, false, ImVec2(200+26*fontSize,300))) {
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
          draw_this_frame = true;

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

      if (ImGui::Button("ImGui Samples")) show_demo_window ^= 1;
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
    compute_ortho_proj_mat(window, rparams.vcx, rparams.vcy, &rparams.vsize, gl_projection);
    sim.drawGL(gl_projection, rparams);

    // if simulation has not been initted, draw the features instead!
    //for (auto const& bf : bfeatures) { bf.drawGL(gl_projection, rparams); }

    // here is where we write the buffer to a file
    if ((is_ready and record_all_frames and sim_is_running) or draw_this_frame) {
      static int frameno = 0;
      std::stringstream pngfn;
      pngfn << "img_" << std::setfill('0') << std::setw(5) << frameno << ".png";
      png_out_file = pngfn.str();
      (void) saveFramePNG(png_out_file);
      std::cout << "Wrote screenshot to " << png_out_file << std::endl;
      frameno++;
      draw_this_frame = false;
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
    // int display_w;
    // int display_h;
    // glfwMakeContextCurrent(window);
    // glfwGetFrameBufferSize(window, &display_w, &display_h);
    // glViewport(0, 0, display_w, display_h);
    // glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    // glClear(GL_COLOR_BUFFER_BIT);
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

