// two versions of a simple cube

#pragma once

#include <cstdint>
#include <vector>

// two tris per side
const std::vector<float> cube0 = {{
0.0, 0.0, 0.0,
0.0, 0.0, 1.0,
0.0, 1.0, 0.0,
0.0, 1.0, 1.0,
1.0, 0.0, 0.0,
1.0, 0.0, 1.0,
1.0, 1.0, 0.0,
1.0, 1.0, 1.0
}};

const std::vector<uint32_t> cube0idx = {{
0,4,5,
0,5,1,
4,6,7,
4,7,5,
1,5,7,
1,7,3,
6,2,3,
6,3,7,
2,0,1,
2,1,3,
2,6,4,
2,4,0
}};


// 8 tris per side - corners all have 6 tris meet there
// THIS IS INCOMPLETE
const std::vector<float> cube1 = {{
0.0, 0.0, 0.0,
0.0, 0.0, 1.0,
0.0, 1.0, 0.0,
0.0, 1.0, 1.0,
1.0, 0.0, 0.0,
1.0, 0.0, 1.0,
1.0, 1.0, 0.0,
1.0, 1.0, 1.0
}};

const std::vector<uint32_t> cube1idx = {{
0,4,5,
0,5,1,
4,6,7,
4,7,5,
1,5,7,
1,7,3,
6,2,3,
6,3,7,
2,0,1,
2,1,3,
2,6,4,
2,4,0
}};

