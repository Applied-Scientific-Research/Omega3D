R"(
#version 150

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
uniform vec4 def_color;
uniform bool use_def;
in vec3 pos;
in float str;
out vec4 base_color;

void main() {
  if (use_def) {
    base_color = def_color;
  } else {
  // flip between base colors based on magnitude of strength
    //base_color = abs(str)*(pos_color*step(0.0f, str) + neg_color*step(0.0f, -str));
    base_color = abs(str)*pos_color;
  }
 
  // pass 1 vert as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(pos, 1.f);
}
)"
