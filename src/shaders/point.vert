R"(
#version 150

uniform mat4 ModelView;
uniform mat4 Projection;
uniform vec4 def_color;
uniform float rad;
in vec4 quad_attr;
in float px;
in float py;
in float posz;
out vec4 base_color;
out vec2 txcoord;

void main() {
  // color pass-through
  base_color = def_color;

  // pass through texture coordinates
  txcoord = quad_attr.xy;

  // make 4 verts as a single primitive and set texture coords - see other shaders
  //gl_Position = Projection * vec4(px + rad*quad_attr.x, py + rad*quad_attr.y, posz, 1.f);

  // do it again, but align the quads to face the camera
  vec4 vpos = vec4(px, py, posz, 1.f);
  gl_Position = Projection * vec4((ModelView*vpos).xyz + rad*quad_attr.xyz, 1.f);
}
)"
