/*
 * Omega3D - The Vorticity Flow Solver
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 * Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Omega3D.h"
#include "Simulation.h"
#include "FlowFeature.h"

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw_gl3.h"
#ifdef _WIN32
  // for glad
  #define APIENTRY __stdcall
  // for C++11 stuff
  #include <ciso646>
#endif
#include "glad.h"
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
    (*_size) *= pow(1.05f, io.MouseWheel);

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
void compute_projection_matrix(GLFWwindow*         _thiswin,
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


// execution starts here

int main(int argc, char const *argv[]) {
  std::cout << std::endl << "Omega3D GUI" << std::endl;

  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  //std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;

  static bool sim_is_running = false;
  static bool begin_single_step = false;
  //const solution_t solver = direct_cpu;
  //std::array<double,Dimensions> fs = {0.0, 0.0, 0.0};
  //double time = 0.0;
  //const double dt = 0.01;

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

  if (!gladLoadGL()) {
    printf("Something went wrong!\n");
    exit(-1);
  }

  // Setup ImGui binding
  ImGui_ImplGlfwGL3_Init(window, true);

  //glfwSetKeyCallback(keyboard_callback);

  // Load Fonts
  // (there is a default font, this is only if you want to change it. see extra_fonts/README.txt for more details)

  // Get and set some IO functions
  ImGuiIO& io = ImGui::GetIO();
  io.IniFilename = ".omega3d.ini";

  // projection matrix for the particles
  float vcx = -0.5f;
  float vcy = 0.0f;
  float vsize = 2.0f;
  std::vector<float> gl_projection;
  compute_projection_matrix(window, vcx, vcy, &vsize, gl_projection);

  // GUI and drawing parameters
  bool show_stats_window = false;
  bool show_terminal_window = false;
  bool show_test_window = false;
  ImVec4 clear_color = ImColor(15, 15, 15);
  ImVec4 pos_circ_color = ImColor(207, 47, 47);
  ImVec4 neg_circ_color = ImColor(63, 63, 255);
  //static bool show_origin = true;
  static bool is_viscous = false;


  // for starters, generate some vortons, particles, and field points
  // eventually use the GUI to generate these elements
  ffeatures.emplace_back(std::make_unique<BlockOfRandom>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.01, 5000));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(10000, active, lagrangian));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(5000, inert, lagrangian));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(2000, inert, fixed));
 
  // Main loop
  while (!glfwWindowShouldClose(window))
  {
    glfwPollEvents();
    ImGui_ImplGlfwGL3_NewFrame();

    //
    // Update simulation
    //

    // get results of latest step, if it just completed
    bool is_ready = sim.test_for_new_results();

    // see if we should start a new step
    if (is_ready and (sim_is_running || begin_single_step)) {

      // if particles are not yet created, make them
      if (not sim.is_initialized()) {
        std::cout << std::endl << "Initializing simulation" << std::endl;

        // initialize particle distributions
        for (auto const& ff: ffeatures) {
          sim.add_particles( ff->init_particles(sim.get_ips()) );
        }

        // initialize solid objects
        //for (auto const& bf : bfeatures) {
        //  sim.add_boundary( bf->get_type(), bf->get_params() );
        //}

        // initialize panels
        //sim.init_bcs();

        sim.set_initialized();
      }

      // generate new particles from emitters
      for (auto const& ff : ffeatures) {
        sim.add_particles( ff->step_particles(sim.get_ips()) );
      }

      // begin a dynamic step: convection and diffusion
      // eventually use the async call here
      sim.step();
      sim.async_step();

      begin_single_step = false;
    }

    // check mouse for drag and rescaling!
    if (not io.WantCaptureMouse) {
      mouse_callback(window, io, &vcx, &vcy, &vsize);
    }

    //
    // The main Omega2D window
    //
    {

      ImGui::Begin("Omega3D");
      ImGui::TextWrapped("Welcome to Omega3D. First set your simulation globals, then add one or more flow structures, then press RUN. Have fun!");

      ImGui::Separator();
      ImGui::Text("Simulation globals");
      ImGui::SliderFloat("Time step", sim.addr_dt(), 0.001f, 1.0f, "%.4f", 1.0f);
      ImGui::Checkbox("Fluid is viscous (diffuses)", &is_viscous);
      if (is_viscous) {
        // show the toggle for AMR
        //static bool use_amr = false;
        //ImGui::Checkbox("Allow adaptive resolution", &use_amr);
        //sim.set_amr(use_amr);
        // and let user choose Reynolds number
        ImGui::SliderFloat("Reynolds number", sim.addr_re(), 1.0f, 10000.0f, "%.1f", 2.0f);
        ImGui::Text("Particle spacing %g", sim.get_ips());
      } else {
        static float my_ips = 0.03141;
        ImGui::SliderFloat("Particle spacing", &my_ips, 0.001f, 1.0f, "%.3f", 2.0f);
        // change underlying Re when this changes
        sim.set_re_for_ips(my_ips);
      }
      ImGui::InputFloat3("Freestream speed", sim.addr_fs());

      ImGui::Separator();
      ImGui::Text("Flow structures");

      ImGui::Separator();
      ImGui::Text("Rendering parameters");
      ImGui::ColorEdit3("positive circulation", (float*)&pos_circ_color);
      ImGui::ColorEdit3("negative circulation", (float*)&neg_circ_color);
      ImGui::ColorEdit3("background color", (float*)&clear_color);
      //ImGui::Checkbox("show origin", &show_origin);

      if (ImGui::Button("Recenter")) {
        // put everything back to center
        vcx = -0.5f;
        vcy = 0.0f;
        vsize = 2.0f;
      }

      ImGui::Separator();
      if (sim_is_running) {
        ImGui::Text("Simulation is running...time = %g", sim.get_time());
        if (ImGui::Button("PAUSE", ImVec2(200,0))) sim_is_running = false;
      } else {
        ImGui::Text("Simulation is not running, time = %g", sim.get_time());
        if (ImGui::Button("RUN", ImVec2(200,0))) sim_is_running = true;
        ImGui::SameLine();
        if (ImGui::Button("Step", ImVec2(120,0))) begin_single_step = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Reset", ImVec2(120,0))) {
        // remove all particles and reset timer
        sim.reset();
      }

      ImGui::Separator();
      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
      //ImGui::Text("Number of panels: %ld  Number of particles: %ld", sim.get_npanels(), sim.get_nparts());
      ImGui::Text("Number of particles: %ld", sim.get_nparts());

      ImGui::End();
    }

    // Rendering
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw the simulation: panels and particles
    compute_projection_matrix(window, vcx, vcy, &vsize, gl_projection);
    sim.drawGL(gl_projection, &pos_circ_color.x, &neg_circ_color.x);

    // draw the GUI
    ImGui::Render();

    // all done! swap buffers to the user can see
    glfwSwapBuffers(window);
  }

  // Cleanup
  ImGui_ImplGlfwGL3_Shutdown();
  glfwTerminate();

  return 0;
}

