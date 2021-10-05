R"(
#version 150
//
// particle.vert - scale particle size and set color
//
// (c)2019 Applied Scientific Research, Inc.
//         Mark J. Stock <markjstock@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

uniform mat4 ModelView;
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
  vec3 str = mat3(ModelView) * vec3(sx, sy, sz);

  // flip between base colors based on whether strength points in or out of the screen
  base_color = pos_color*step(0.0f, str.z) + neg_color*step(0.0f, -str.z);
  //base_color = pos_color*step(0.0f, sz) + neg_color*step(0.0f, -sz);
  //base_color = vec4(1.f, 1.f, 1.f, 1.f);

  // pass through texture coordinates
  txcoord = quad_attr.xy;

  // magnitude of strength density scales fragment color
  //strength = sqrt(sx*sx + sy*sy + sz*sz) / (r*r*r);

  // or magnitude of strength normal to screen?
  // note: have to keep this next line or else glsl will compile out the sx and sy vars and crash!
  //strength = sx + sy + sz;
  //strength = abs(sz) / (r*r*r);
  strength = abs(str.z) / (r*r*r);

  // make 4 verts as a single primitive and set texture coords - see other shaders
  float rscale = 2.5f * r * rad_scale;
  //gl_Position = Projection * vec4(px + rscale*quad_attr.x, py + rscale*quad_attr.y, posz, 1.f);

  // do it again, but align the quads to face the camera
  vec4 vpos = vec4(px, py, posz, 1.f);
  gl_Position = Projection * vec4((ModelView*vpos).xyz + rscale*quad_attr.xyz, 1.f);
}
)"
