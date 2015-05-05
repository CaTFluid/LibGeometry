// LibMaterial,  created 04/04/2014 Wenjun Kou
// see example in libmesh/system_of_equations/example 6


// 1) info of materialial parts: subdomian id, (material id)

#ifndef INC_LIBMATERIAL_ENUM
#define INC_LIBMATERIAL_ENUM
#include <string>

namespace LibMaterial{
static const double PI = 3.1415926;
static const std::string ErrorMsg= "[Error:]";
static const std::string WarnMsg= "[Warning:]";
  
  
  
enum  LayerID{
   Layer_1= 0, // inner mucosal layer
   Layer_2, Layer_3, Layer_4, Layer_5, Layer_6, Layer_7
  
}; // LayerID
  
enum TubeBoundaryID{ // for tubular structure
 BOUNDARY_ID_top = 0,
 BOUNDARY_ID_bot,
 BOUNDARY_ID_inner,
 BOUNDARY_ID_outer
 
  
}; // BoundaryID
enum  BoundaryTypeID {
  BOUNDARY_TYPE_Dirichlet =0,
  BOUNDARY_TYPE_TractionForce, // we specify the F1, F2, F3
  BOUNDARY_TYPE_PRESSURE, // we speficy p, but need to get n and then p1,p2,p3
  BOUNDARY_TYPE_FRICTION // we specify friction: two f1, f2, need to compute the normal
  

};
}

using namespace LibMaterial;
 
#endif // LibMaterial