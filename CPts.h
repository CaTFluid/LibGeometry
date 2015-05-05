// class of point operation, not aim for array of storage
// current 3d case,
#ifndef _INCL_CPTS
#define _INCL_CPTS
#include <iostream>
#include <cmath>
using namespace std;
//#include <cmath>
namespace MY_APP
{ 
  class CPts{
   /// PUBLIC /// mainly for read/write
public:
  /// cint >> pt (notice operator is outside-call)
  // thus to access private members, we add friend
  
  friend istream& operator>> (istream& is, CPts& pt)
  { is>>pt.x>>pt.y>>pt.z;
   return is;
  }
  /// <<
  friend ostream& operator<<(ostream& os, const CPts& pt)
  { os<<pt.x<<"\t"<<pt.y<<"\t"<<pt.z;
    return os;
  }
  /// +  cpt1+ cpt2
  CPts operator+(const CPts &rhs)
  { return CPts(x+rhs.x,y+rhs.y,z+rhs.z);
  }
  
  /// - cpt1-cpt2
  CPts operator-(const CPts &rhs)
  {return CPts(x-rhs.x,y-rhs.y,z-rhs.z);
  }
  ///// constructor 
  CPts()
  {x=0.0;y=0.0;z=0.0;
  }
  // given cods
  CPts(double dx, double dy, double dz)
    :x(dx),
     y(dy),
     z(dz)
     {}
     // not reference, just copy the value
   CPts(const CPts& that)
   { x=that.x;
     y=that.y;
     z=that.z;
   }
   //// ulility:
   // distance

  double x;
  double y;
  double z;
};
//// outside call + ,- 
inline
CPts operator-(const CPts &left, const CPts &right)
{ return CPts(left.x-right.x, left.y-right.y, left.z-right.z);
}
inline
CPts operator+(const CPts &left, const CPts &right)
{ return CPts(left.x+right.x, left.y+right.y, left.z+right.z);
}
//// utility function for points
//1) distance of 2 points
inline   
double dist(CPts pt1, CPts pt2)
   { double dx=pt1.x-pt2.x;
     double dy=pt1.y-pt2.y;
     double dz=pt1.z-pt2.z;
     return sqrt(dx*dx+dy*dy+dz*dz);
   }// distance 
// 2) dot product of 2 vector(points)
inline
double dot_prod(CPts vec1, CPts vec2)
{return (vec1.x*vec2.x+vec1.y*vec2.y+vec1.z+vec2.z);
}// dot product
// 3) cross_prod of 2 vector(points) return a vector
inline
CPts cross_prod(CPts vec1, CPts vec2)
{double my_x=vec1.y * vec2.z-vec2.y*vec1.z; //y1*z2-y2*z1
 double my_y=vec2.x * vec1.z -vec1.x * vec2.z;// x2*z1 - x1*z2
 double my_z=vec1.x * vec2.y -vec2.x * vec1.y;// x1*y2 - x2*y1
 return CPts(my_x,my_y,my_z);
}
// 4) tri_prod: z1 x z2 . z3
inline 
double tri_prod(CPts vec1, CPts vec2, CPts vec3)
{ return dot_prod(vec3,cross_prod(vec1,vec2));
}
// 5) norm: norm_2 of a vector
inline
double norm_2(CPts vec)
{return sqrt(dot_prod(vec,vec));
}
   
   
}


#endif
