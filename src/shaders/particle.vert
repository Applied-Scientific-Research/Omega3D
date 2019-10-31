R"(
#version 150

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
uniform float rad_scale;
in vec4 quad_attr;
in float px;
in float py;
in float posz;
in float sx;
in float sy;
in float sz;
in float r;
out vec4 base_color;
out vec2 txcoord;
out float strength;

void main() {

  // find strength vector in projection
  //vec4 str = Projection * vec4(sx, sy, sz, 1.f);

  // flip between base colors based on whether strength points in or out of the screen
  //base_color = pos_color*step(0.0f, str.z) + neg_color*step(0.0f, -str.z);
  base_color = pos_color*step(0.0f, sz) + neg_color*step(0.0f, -sz);
  //base_color = vec4(1.f, 1.f, 1.f, 1.f);

  // pass through texture coordinates
  txcoord = quad_attr.xy;

  // magnitude of strength density scales fragment color
  //strength = sqrt(sx*sx + sy*sy + sz*sz) / (r*r*r);

  // or magnitude of strength normal to screen?
  // note: have to keep this next line or else glsl will compile out the sx and sy vars and crash!
  strength = sx + sy + sz;
  strength = abs(sz) / (r*r*r);

  // make 4 verts as a single primitive and set texture coords - see other shaders
  float rscale = 2.5f * r * rad_scale;
  gl_Position = Projection * vec4(px + rscale*quad_attr.x, py + rscale*quad_attr.y, 0.f, 1.f);
}
)"
