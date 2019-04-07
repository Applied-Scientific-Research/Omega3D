R"(
#version 150

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
uniform vec4 def_color;
in float px;
in float py;
in float posz;
in float sx;
in float sy;
in float sz;
out vec4 base_color;
out float strength;

void main() {
  // flip between base colors based on magnitude of strength
  //base_color = pos_color*step(0.0f, rawstr) + neg_color*step(0.0f, -rawstr);
  base_color = def_color;

  // magnitude of strength scales fragments
  //strength = 1.f;
  strength = sqrt(sx*sx+sy*sy+sz*sz);

  // make 4 verts as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(px, py, posz, 1.f);
}
)"
