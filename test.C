// test.C to
// conclusion: better to put all the definition together, usin inline function
// operator >> over load for Class and Structor
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
// class CPts
class CPts{
   /// PUBLIC /// mainly for read/write
public:
  /// cint >> pt (notice operator is outside-call)
  // thus to access private members, we add friend
  double x;
  double y;
  double z;
  friend istream& operator>> (istream& is, CPts& pt)
  { is>>pt.x>>pt.y>>pt.z;
   return is;
  }
  /// <<
  friend ostream& operator<<(ostream& os, const CPts& pt)
  { os<<pt.x<<"\t"<<pt.y<<"\t"<<pt.z;
    return os;
  }
};
// structure Cnodes
struct SPts
{ double x;
  double y;
  double z;

};
struct vecSpts
{ int num;
  vector <SPts> vec_pts;
};

// output
ostream& operator<<(ostream& os, const SPts& pt)
 { os<<pt.x<<"\t"<<pt.y<<"\t"<<pt.z;
    return os;

 }

ostream& operator<<(ostream& os, const  vecSpts& vecpt)
 { int num = vecpt.num;
   for (int k=0; k < num; k++)
   { os << vecpt.vec_pts[k];   
    }
	return os;
 }

int main ()
{ double a, b, c;
  a= 1.0; b= 0.5*a; c= 0.2*a;
  SPts S_pt;
  S_pt.x =a;
  S_pt.y =b;
  S_pt.z =c;
  CPts C_pt;
  C_pt.x =a;
  C_pt.y =b;
  C_pt.z =c;
 cout << "print Struc pts: " << S_pt << endl;
 cout << "print Class pts: " << C_pt << endl;
 // try vector
 vector<SPts> vec_spt;
 vec_spt.resize(1);
 vec_spt[0].x =a;
 vec_spt[0].y =b;
 vec_spt[0].z =c;
 cout << "print vector pts: " << vec_spt[0] << endl;
// try vector structure
 vecSpts vec_newspt;
 vec_newspt.num = 2;
 vec_newspt.vec_pts.resize(2);
 vec_newspt.vec_pts[0].x =a;
 vec_newspt.vec_pts[0].y =b;
 vec_newspt.vec_pts[0].z =c;
 cout << "print new sturcture vector pts: " <<  vec_newspt<< endl;
}
