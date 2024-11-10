#include "Unstructured_FV_Solver.h"
#include "Unstructured_Mesh.h"


///////////////////////////////////////////////////////////////////////
//                         Main
///////////////////////////////////////////////////////////////////////
int main() {
 
    auto solver = Unstructured_FV_Solver("nofins.msh");
    solver.make_movie(0.2,0.002,401,"Forced_no_fins_movie");

  return 0;
}


