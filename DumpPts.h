// concrete class for read vtk file and dump points to related file
// created 04/05/2013, refer to the TimeDependentSpring
// Operator members include
// 1) read input context, 2) read vtk, dump it according to the context as back-up data ( level by level)
// 3) calculate area/.. and dump the result
#ifndef _INCL_DUMP_PTS
#define _INCL_DUMP_PTS
// C++ STDLIB INCLUDES
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>
#include <map>

#include "CPts.h"
namespace MY_APP
{
static const double PI = 3.1415926535897932384626433832795028;
// function-ptr type for calculate volume at each file(time snap)
// return the volume, given the array of area at different levels
typedef // vec_para: area, z_level and so on;
void (*VolumeFucPtr)(const std::vector<std::vector<double> >& vec_para 
		     ,const int Num_levels
		     ,double& return_volume);
// default: A1-A2 (dz)= dv=
inline void
default_volume_fcn(const std::vector<std::vector<double> >& vec_para 
		     ,const int Num_levels
		     ,double& return_volume)
{
  return_volume=0.0;
  int it_level;
  for(it_level=0;it_level<Num_levels-1;it_level++)
  { double A1=vec_para[it_level][0];
    double A2=vec_para[it_level+1][0];
    double dz=vec_para[it_level+1][1]-vec_para[it_level][1];
    if (dz<0.0) dz=0.0-dz;
    return_volume+=1.0/3.0*(A1+A2+sqrt(A1*A2))*dz;
    
  }
}








// function-ptr type for calculate area on each level, based on coordinates
typedef // return the area, also we could modify the parameters for dump
void (*AreaFucPtr)(const std::vector<CPts>& Level_points, const int Num_points, 
		     std::vector<double>& para, int& Num_para);


// Amneet's method: default fun-ptr to calculate area_fcn_ptr
// return: area z-cod, min-d, max_d, minnode_id, maxnode_id;
inline void 
default_area_fcn(const std::vector<CPts>& Level_points, const int Num_points,
		  std::vector<double>& para, int & Num_para)
{ //1) calulate the center of level, assume x=y=0;, we only need to compute z;
   // as we need to compute max, min-distance, and level cods
  // get the level and center of the area: average z coordinate
  int it_pt;
  double level_z=0.0;
  for(it_pt=0; it_pt<Num_points; it_pt++)
     level_z+=Level_points[it_pt].z;
  level_z=level_z/Num_points; // get averge 
  CPts	center(0,0,level_z);
  int max_id;
  int min_id;
  CPts max_pt;
  CPts min_pt;
  double max_dist=0.0;
  double min_dist=dist(Level_points[0],center);
  //
  for (it_pt=0;it_pt<Num_points;it_pt++)
  {
    double cur_dist=dist(Level_points[it_pt],center);
    if(cur_dist>max_dist)
    { max_dist=cur_dist;
      max_id=it_pt;
      max_pt=Level_points[it_pt];
    }
    if (cur_dist<min_dist)
    {
      min_dist=cur_dist;
      min_id=it_pt;
      min_pt=Level_points[it_pt];
    }
  }
  //
  double area=PI*min_dist*max_dist;
  Num_para=6;
  para.resize(Num_para);
  para[0]=area;
  para[1]=level_z;
  para[2]=min_dist;
  para[3]=max_dist;
  para[4]=min_id;
  para[5]=max_id;
  return;
  
}

// begin the definition of the class
class DumpPts
{
public:
    // Constructor
    DumpPts( std::string InfoFile);
    
    
    // Dump pts to file by the order of level (top 2 Bot) 
    // input: file_id, by which we know which file we read
    //  
    bool
    DumpPtsBot2Top(int file_id);
    bool
    DumpAreaBot2Top(int file_id, vector<double> & volume);   
    
    
    // version 4: 09/05/2013-- dump spring rest length
    // this will be called in DumpPtsBot2Top
    bool
    DumpLengthBot2Top(const std::vector<CPts> & vec_Pts, int file_id);
    
    void
    RegisterComputAreaFcn(int fcn_idx, const AreaFucPtr area_fcn_ptr, const std::string fcn_info
			  , const VolumeFucPtr volume_fcn_ptr)
    {d_area_fcn_map[fcn_idx]=area_fcn_ptr;
     d_fcn_info_map[fcn_idx]=fcn_info; // to give information on the ouput of fcn;
     d_volume_fcn_map[fcn_idx]=volume_fcn_ptr;
    }
    // get total number of file we need to deal with
    int GetTotalNum(){return d_total_files;}
    std::string GetOuputPrefix() {return d_oup_prefix;}
    int GetLayerNum(){return d_num_layers;}
    
private:
     //>> mesh_info
    int d_num_axial; // number of axial levels
    int d_num_circ; // number of circular nodes
    //>> node_id inof
    int d_id_bgn;
    int d_id_end;
    int d_fcn_id;
    //>> dump Pts format:
    // number of points per row; must be dividable by the num of d_num_circ;
    int d_num_pts_per_row;
    //>> file_info
    std::string d_inp_prefix; // prefix of input filename
    std::string d_oup_prefix; // prefix of output filename
    int d_total_files; // total number of files
    // we provide interface for area function 
    std::map<int, AreaFucPtr> d_area_fcn_map;
    std::map<int, std::string> d_fcn_info_map;
    std::map<int, VolumeFucPtr> d_volume_fcn_map;
    //>> add function of dumping multiple layers, 05/23/2013
    int d_version; // version_number
    int d_num_layers;// number of layers
    std::vector<int> d_file_ids; // non_continus id;
    //<< 05/23/2013
    //>> 09/05/2013
    int d_is_dump_rad_spring;
    int d_is_dump_circ_spirng;
    int d_is_dump_axi_spring;
    //<< 09/05/2013
    
};// class DumpPts


}// MY_APP
#endif