// version 5:   12/17/2013
//		compatiable with previous case: dump all files, need to use loop
//		current case, we only deal with continous files: t000 t001... t010
//		line 5: use "-1" as marker of continous file input
//		-1	start_id(default =0)	end_id = start_id + file number -1;	
// version 4:   09/05/2013
//		dump points --> spring length, 
//		dump spring length --> spring force, according to parameters 

// version 3: dump multiples layers, with non-continus file_id
//	05/23/2013 (IM_t*_layer*.pts; .area. volume
// concrete class for read vtk file and dump points to related file
// version 2: the start node_id may begin from non-zero, 04/15/2013 >>
//		we may need to deal with multiple structures, <<
// version 1:created 04/05/2013, refer to the TimeDependentSpring
// Operator members include
// 1) read input context, 2) read vtk, dump it according to the context as back-up data ( level by level)
// 3) calculate area/.. and dump the result

// C++ STDLIB INCLUDES
/*
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <math.h>
#include <map>
*/
#include "DiscardComments.h"
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
// include pts class
#include "DumpPts.h"
//#include "CPts.h"

using namespace MY_APP;
using namespace MY_CPP;
// constructor
DumpPts::DumpPts(std::string InfoFile)
{
  std::string line_string;
  std::ifstream file_stream;
  std::string FileName=InfoFile;
  std::string OutPrefix="[MESH_INFO]:";
  file_stream.open(FileName.c_str(),std::ios::in);//input(to the screen)=read
  if (file_stream.is_open())
 { 
    int it_line;
    int total_lines=6; // added 05/23/2013 for multiples layers; line 5: list of file_id
  for (it_line=0; it_line< total_lines;it_line++)
  {
      if(!std::getline(file_stream,line_string))
	// error
      {std::cout<< "error for read"<<std::endl;
	return ;
	
      }
      else
	{  
	  
	  // first output the file info
	  std::cout<< OutPrefix << line_string << std::endl;
	   line_string=discard_comments(line_string);          
           std::istringstream line_stream(line_string);
	  if (it_line==0){
	   //1) 1st line: start_id, end_id, and fcn_ptr
	  line_stream>>d_id_bgn >> d_id_end >>d_fcn_id>>d_version;
	  std::cout << OutPrefix <<"Check:id-bgn, id_end, fcnid: "
	    << d_id_bgn <<", "<< d_id_end<<", "<<d_fcn_id<<std::endl;
	  std::cout <<OutPrefix <<"!!Notice: we use version: " << d_version << std::endl;
	  if (d_id_end < 1)
	  {
	    if (d_version<3){
		std::cout<< OutPrefix <<"Warn: id_end(default)= axial_no * num_circ+ id_bgn-1"
			<< std::endl;}
	    else{std::cout<< OutPrefix <<"Warn: id_end(default)= axial_no * num_circ*num_layers+ id_bgn-1"
			<< std::endl;}
	    }
	      
	
	  }
	  if (it_line==1)
	  { //2) 2nd line; mesh_info, num_axial, num_circ
	  line_stream >>d_num_axial >> d_num_circ >> d_num_layers;
	   std::cout << OutPrefix <<"Check:num_aixal, num_circ, num_layers: "
	    << d_num_axial <<", "<< d_num_circ << "," << d_num_layers<< std::endl;
	  if (d_id_end < 1){
	    if (d_version<3)
		{ d_id_end=d_id_bgn+d_num_axial * d_num_circ -1;
		std::cout<< OutPrefix <<"Check:id_end: "<< d_id_end << std::endl;
		}
	    else {d_id_end=d_id_bgn+d_num_axial * d_num_circ * d_num_layers -1;
		std::cout<< OutPrefix <<"Check:id_end: "<< d_id_end << std::endl;
	      
	    }
		  
	  }
	  }
	  if (it_line==2)
	  {//3) prefix for the file: 
	  line_stream >> d_inp_prefix >> d_oup_prefix;
	  std::cout << OutPrefix <<"Check:inp_file_prefix, oup_file_prefix: "
	    << d_inp_prefix <<"***.vtk; "<< d_oup_prefix 
	    << "***.pts,***.area"<<std::endl;
	  }
	  if (it_line==3)
	  { //4) output info: file numbers and num.pts per row in each file
	  line_stream >> d_total_files >> d_num_pts_per_row;
	  std::cout << OutPrefix <<"Check:output_fileNums, PtsPerRow: "
	    << d_total_files <<", "<< d_num_pts_per_row <<std::endl;
	  }
	  if (it_line==4) // added for version 3: 05/23/2013
	  { //5) files id
	    d_file_ids.resize(d_total_files);
	    int fid=0;
	    std::cout<<OutPrefix <<"Check files:" << std::endl;
	// changed for version 5: 12/17/2013 >>
	    int IdMarker;
	    line_stream >> IdMarker;
	  	
            if (IdMarker < 0) { // continuous file id
		std::cout<<OutPrefix <<"Continuous files, starting :";
		line_stream >> d_file_ids[0];
		std::cout << d_file_ids[0];
		for (fid=1;fid<d_total_files;fid++)
	    		{ 
			  d_file_ids[fid]=d_file_ids[fid-1]+1;
	    		}
		
	      }
	   else { // several files supposed
		std::cout<<OutPrefix <<"Given several files:";
		d_file_ids[0]=IdMarker;
		std::cout << d_file_ids[0] << ", ";
		if (d_total_files > 1) 
		{
		   for (fid=1;fid<d_total_files;fid++)
	    	     {	line_stream >> d_file_ids[fid];
	      		std::cout << d_file_ids[fid] << ", ";
	    	     }
		}
	    }	
	/*    for (fid=0;fid<d_total_files;fid++)
	    {line_stream >> d_file_ids[fid];
	      std::cout << d_file_ids[fid] << ", ";
	    }
	 */
	// changed for version 5: 12/17/2013 <<
	    std::cout << std::endl;
	  }
	  if (it_line == 5) // added for version 4: 09/05/2013
	  { // bool id for dump springs
	    d_is_dump_rad_spring=0;
	    line_stream >> d_is_dump_rad_spring;
	    std::cout<< OutPrefix << "Dump radial spring length?" << ((d_is_dump_rad_spring>0) ? "Yes" : "no") <<std::endl;
	    d_is_dump_circ_spirng=0;
	    line_stream >> d_is_dump_circ_spirng;
	    std::cout<< OutPrefix << "Dump circ spring length?" << ((d_is_dump_circ_spirng>0) ? "Yes" : "no" )<<std::endl;	    
	    d_is_dump_axi_spring=0;
	    line_stream >> d_is_dump_axi_spring;
	    std::cout<< OutPrefix << "Dump axial spring length?" << ((d_is_dump_axi_spring>0) ? "Yes" : "no") <<std::endl;		    
	  }
	}
    
  }
      file_stream.close();
 }
   // register default area_fcn_ptr
   std::string default_fcn_info="area;z_level;min_dist;max_dist;min_nodeId;max_nodeId";
   RegisterComputAreaFcn(0, &default_area_fcn, default_fcn_info,&default_volume_fcn); 
}// Constructor

// given file_id, we read the pts and write it into file
bool DumpPts::DumpPtsBot2Top(int file_id)
{ // check file_id, and get output filename
  bool IsReadDone=false;
  bool IsWriteDone=false;
  int total_pts=d_id_end-d_id_bgn+1;
  std::vector<CPts> vector_pts(total_pts);
  std::string inpFile;
  std::string oupFilePts;
  std::string OutPrefix="[PTS_INFO]:";
  if (file_id>=d_total_files)
  {std::cout<< OutPrefix<< "error on file-id, too big" << std::endl;
   return false;
  }
  else
  { // 0) version 3, discontunuous file_id >> 05/23/2013
    if(d_version>=3)
    { if (d_file_ids.size()-1<file_id) // error
      {std::cout<< OutPrefix<< "error: ver 3:  file_id > size of d_files_ids" << std::endl;
      }
      else // get_new file_id based on the file_ids
      {file_id=d_file_ids[file_id];
      }
    } // << 05/23/2013
    // 1) get the related filename:
    // append middle part: ***
    char temp[10];
    sprintf(temp,"%03d", file_id);
    inpFile=d_inp_prefix+temp+".vtk";
    
    
    // 2) read file
  std::string line_string;
  std::ifstream file_stream;
  //std::string FileName=InfoFile;
  file_stream.open(inpFile.c_str(),std::ios::in);
  if(file_stream.is_open())
  { 
   int it_pt=0;
      // 2.1) cout comment lines ;WE KNOW first 8 lines are commented
   int num_comment_lines=8;
   int it_line;
   std::cout<< "+++++++++++++++++++++++++++++++below is vtk_header-info" << std::endl;
   for (it_line=0;it_line<num_comment_lines;it_line++)
   { if(!std::getline(file_stream,line_string))
      {std::cout << "error for read line:" << it_line << std::endl;
	return false;
      }
      else{
	std::cout << line_string << std::endl;
      }
   }
   std::cout<< "+++++++++++++++++++++++++++++++below is to get points num=" << total_pts << std::endl;
    // 2.2) read data to the big array
   // version 2: deal with non-zero d_start_id 04/15/2013 >>
   // ignore nodes before the start id
   int id_ignore=0;
   CPts ignore_pt;
   while (id_ignore<d_id_bgn)
   {
     file_stream >> ignore_pt;
     id_ignore++;
   }
   // 04/15/2013 <<
   
   for (it_pt=0; it_pt<total_pts; it_pt++)
   {	
     file_stream >> vector_pts[it_pt];
   }
    IsReadDone=true;
    file_stream.close();
  }
  else
  {std::cout<<"can not open the file"<<std::endl;
   return false;
  }
  }
  //3) output/write pts
  //// 3.1) prepare infomation:
  if (!IsReadDone)
    return false;
  else
    std::cout << OutPrefix<< "Done: Read vtk file:" << inpFile<<std::endl;
  //////////////////////////////////
  //>> added 05/23/2013, multiple layers
  std::string layer_fix="_layer";
  char temp[10];
  sprintf(temp,"%03d", file_id);
  int ind_layer;
  for (ind_layer=0;ind_layer<d_num_layers; ind_layer++)
  {
    char strlayer[10];
    sprintf(strlayer,"%01d", ind_layer);
  oupFilePts=d_oup_prefix+temp+layer_fix+strlayer+".pts";
  std::ofstream file_stream;
  file_stream.open(oupFilePts.c_str(),std::ios::out);
  std::cout << OutPrefix<<"Begin: Write pts file:" << oupFilePts << std::endl;
  //// 3.2) write general information
  //>> level_no, circ_no, and num_pts per row
  file_stream << d_num_axial << "\t" << d_num_circ 
    << "\t" << d_num_pts_per_row << "\t" << ind_layer
    << "\t #num_aixal level, num_circ pts, num_pts per row, ind_layer"
    << "\n";
  int layer_offset=ind_layer*d_num_axial*d_num_circ;
  int it_level;
  int it_pt;
  // loop1) each level, top to bottom
  for(it_level=0; it_level<d_num_axial; it_level++ )
  { // and the array-id interval=d_num_axial;
    // >> change 05/23/2013
    //int array_id_start=it_level;
      int array_id_start=it_level+layer_offset;
    // << change 05/23/2013
    for (it_pt=0; it_pt<d_num_circ; it_pt++)
    { file_stream << vector_pts[array_id_start+it_pt*d_num_axial] << "\t";
      if ((it_pt+1) % d_num_pts_per_row ==0 )// add newline
	file_stream << "\n";
      
      
    }
  }
  file_stream.close();
  IsWriteDone=true;
  
  
    
  
  //vector_pts.clear();
  std::cout << OutPrefix<<"Done: write to .pts file from bottom to top:" << oupFilePts<<std::endl;
  
  }// loop of ind_layer
  
  // added 09/05/2013 dump length
  
 if( ! DumpLengthBot2Top(vector_pts,file_id))
 {std::cout << "error on dumping length of springs" << std:: endl;}
  
  return IsWriteDone;


}// Dumpt pts to file order: bottom to top;


bool DumpPts::DumpLengthBot2Top(const vector< CPts >& vec_Pts, int file_id)
{
std::string OutPrefix="[SPRING_INFO]:";
/*
 if (file_id>=d_total_files)
  {std::cout<< OutPrefix<< "error on file-id, too big" << std::endl;
   return false;
  }
  else
  { // 0) version 3, discontunuous file_id >> 05/23/2013
    if(d_version>=3)
    { if (d_file_ids.size()-1<file_id) // error
      {std::cout<< OutPrefix<< "error: ver 3:  file_id > size of d_files_ids" << std::endl;
      return false;
      }
      else // get_new file_id based on the file_ids
      {file_id=d_file_ids[file_id];
      }
    } // < 
  }
  */
 //>> OutPut layer by layer; notice radial spring; outer-layer =0; axial spring:top=0; 
  std::string layer_fix="_layer";
  char temp[10];
  sprintf(temp,"%03d", file_id);
  int ind_layer;
  for (ind_layer=0;ind_layer<d_num_layers; ind_layer++)
  { // deal with file name
    char strlayer[10];
    sprintf(strlayer,"%01d", ind_layer);
    string oupFilePre=d_oup_prefix+temp+layer_fix+strlayer; 
    std::vector< std::string> outFileSprings(3,oupFilePre);
    outFileSprings[0]=outFileSprings[0]+".RSpring";
    outFileSprings[1]=outFileSprings[1]+".CSpring";
    outFileSprings[2]=outFileSprings[2]+".ASpring";
    // open the file
    std:: vector< std::ofstream *> vec_file_stream;
    vec_file_stream.resize(3);
    for (int it_type=0; it_type<3; it_type++)
    {
    vec_file_stream[it_type]=new std::ofstream(outFileSprings[it_type].c_str(),std::ios::out);
    *(vec_file_stream[it_type]) << d_num_axial << "\t" << d_num_circ 
    << "\t" << d_num_pts_per_row << "\t" << ind_layer
    << "\t #num_aixal level, num_circ pts, num_pts per row, ind_layer"
    << "\n";
    }
    // get the start id, for current layer:
    std::cout<< "open the file for output: succeed" << std::endl;
    int layer_offset=ind_layer*d_num_axial*d_num_circ;
    int layer_interval=d_num_axial*d_num_circ;
    int circ_interval=d_num_axial;
    
     // loop1) each level, top to bottom
    int it_level;
    int it_pt;
    for(it_level=0; it_level<d_num_axial; it_level++ )
    { // and the array-id interval=d_num_axial;
    // >> change 05/23/2013
    //int array_id_start=it_level;
      int array_id_start=it_level+layer_offset;
    // << change 05/23/2013
    for (it_pt=0; it_pt<d_num_circ; it_pt++)
    {   // first node:
        int node_ref_id=array_id_start+it_pt*d_num_axial;
	
	// radial spring
	int node_rad_id = (ind_layer < d_num_layers-1) ? (node_ref_id + layer_interval) : node_ref_id;
        double rLen= dist(vec_Pts[node_ref_id],vec_Pts[node_rad_id]);
	if (d_is_dump_rad_spring>0) {
	  *(vec_file_stream[0]) << rLen << "\t";
	  if ((it_pt+1) % d_num_pts_per_row ==0 )// add newline
	  *( vec_file_stream[0]) << "\n";
	}
	
	// circular spring
	int node_circ_id =(it_pt<d_num_circ-1) ? (node_ref_id+circ_interval) : (node_ref_id+circ_interval-layer_interval); 
	double cLen = dist(vec_Pts[node_ref_id], vec_Pts[node_circ_id]);
	if (d_is_dump_circ_spirng>0) {
	 *(vec_file_stream[1]) << cLen << "\t";
	 if ((it_pt+1) % d_num_pts_per_row ==0 )// add newline
	   *(vec_file_stream[1]) << "\n";
	  
	}
	
	// axial spring
	int node_axi_id= (it_level < d_num_axial-1) ? (node_ref_id+1): node_ref_id;
	
	double aLen= dist(vec_Pts[node_ref_id],vec_Pts[node_axi_id]);
	if (d_is_dump_axi_spring) {
	  *(vec_file_stream[2]) << aLen << "\t"; 
	  if ((it_pt+1) % d_num_pts_per_row ==0 )// add newline
	   *(vec_file_stream[2]) << "\n";
	  
	}
        
      

      
      }
    }
    for (int it_type=0; it_type<3; it_type++)
    {vec_file_stream[it_type]->close();
      free(vec_file_stream[it_type]);
    }
    
    std::cout << OutPrefix<<"Done: write to 3 spring file from bottom to top:" << oupFilePre<<std::endl;
  }
  
  return true;
 
}


bool DumpPts::DumpAreaBot2Top(int file_id, vector<double>& volume)
{
  bool IsReadDone=false;
  bool IsWriteDone=false;
 // int total_pts=d_id_end-d_id_bgn+1;
 // std::vector<CPts> vector_pts(total_pts);
  
  
  int total_levels=d_num_axial;
  //std::vector<double> vector_area(total_levels);

  
  
  //// part 1 read data;
  std::string inpFile;
  std::string oupFileArea;
  std::string OutPrefix="[PTS_INFO]:";

   if (file_id>=d_total_files)
  {std::cout<< OutPrefix<< "error on file-id, too big" << std::endl;
   return false;
  }
  // begin to read each layer the write each layer
  //>> version 3, change file_id
  file_id=d_file_ids[file_id];
  
  char temp[10];
  sprintf(temp,"%03d", file_id);
  std::string layer_fix="_layer";
  
  // >> dump summary data for analysis: area, axial-length; (all layers)
  std::vector<std::vector<double> > vec_summary(d_num_axial);
  int ind_axial;
  for(ind_axial=0;ind_axial<d_num_axial;ind_axial++)
  {vec_summary[ind_axial].resize(d_num_layers*2); // area, axial-length
  }
  string oupFileSum=d_oup_prefix+temp+".sum";
  int ind_layer;
  for (ind_layer=0;ind_layer<d_num_layers; ind_layer++)
  {
  
  std::vector< std::vector<double> > vector_para;
  std::vector< std::vector<CPts> > vector_pts;
  vector_para.resize(total_levels);
  vector_pts.resize(total_levels);
    char strlayer[10];
    sprintf(strlayer,"%01d", ind_layer);
  inpFile=d_oup_prefix+temp+layer_fix+strlayer+".pts";
  oupFileArea=d_oup_prefix+temp+layer_fix+strlayer+".area";
  
  std::string line_string;
  std::ifstream file_stream1;
  //std::string FileName=InfoFile;
  file_stream1.open(inpFile.c_str(),std::ios::in);
   if(file_stream1.is_open())
  { if(std::getline(file_stream1, line_string))
    { // just to compare
    std::cout << OutPrefix<< line_string << std::endl;
    std::cout<< OutPrefix << "cur_mesh info:" << d_num_axial <<"\t"
    << d_num_circ << "\t" << d_num_pts_per_row << "\t"<<ind_layer<<std::endl;
    }
    // then we begin to read all data at once
    std::cout<< OutPrefix<< "Begin: read pts in file:" << inpFile <<std::endl;
    int it_level;
    for (it_level=0; it_level<total_levels; it_level++)
    { vector_pts[it_level].resize(d_num_circ);
      int it_pts;
      for (it_pts=0; it_pts<d_num_circ; it_pts++)
	      file_stream1 >> vector_pts[it_level][it_pts];
    }
    // then finish reading
    file_stream1.close();
    IsReadDone=true;
    std::cout<< OutPrefix<<"Done: read pts in file:" << inpFile << std::endl;    
  }
  else 
  { std::cout << OutPrefix<< "can not open file:" << inpFile << std::endl;
    return false;
  }
  
  ///// part 2 calculate and write data
  if (!IsReadDone)
    return false;
  std::cout<< OutPrefix <<"Begin: calculate area with fcn_ptr:"<< d_fcn_id << std::endl;
  int it_level;  
  int num_para;
  for (it_level=0; it_level<total_levels; it_level++)
  {
    (d_area_fcn_map[d_fcn_id])(vector_pts[it_level],d_num_circ
      ,vector_para[it_level],num_para);
 
  }
  ////2.1) finish the calculation
  std::cout<< OutPrefix <<"Done: calculate area with fcn_ptr:"<< d_fcn_id << std::endl;
  ////2.1.2) calculate axial spring-length over each level: average, min, max >> added 05/24/2013
  std::vector<std::vector<double> > vec3_axial_length(total_levels);
  std::string oupFileDisp=d_oup_prefix+temp+layer_fix+strlayer+".disp";
  for (it_level=0;it_level<total_levels-1;it_level++)
  { double sum_length=0.0;
    double min_length=10.0;
    double max_length=0.0;
    int it_pts;
    for (it_pts=0;it_pts<d_num_circ;it_pts++)
    { double cur_length=vector_pts[it_level+1][it_pts].z-vector_pts[it_level][it_pts].z;
      sum_length+=cur_length;
      min_length=cur_length<min_length ? cur_length : min_length;
      max_length=cur_length>max_length ? cur_length : max_length;
    }
     vec3_axial_length[it_level].resize(3);
     vec3_axial_length[it_level][0]=sum_length/d_num_circ;
     vec3_axial_length[it_level][1]=min_length;
     vec3_axial_length[it_level][2]=max_length;
  }
  // the top level: just copy total-levels-1 data
  vec3_axial_length[total_levels-1].resize(3);
  vec3_axial_length[total_levels-1][0]=vec3_axial_length[total_levels-2][0];
  vec3_axial_length[total_levels-1][1]=vec3_axial_length[total_levels-2][1];
  vec3_axial_length[total_levels-1][2]=vec3_axial_length[total_levels-2][2];
  std::cout<< OutPrefix <<"Done: calculate axial displacement" << std::endl;
  std::cout << OutPrefix<<"Begin: write axial disp. in file:" << oupFileDisp<< std::endl;
  std::ofstream file_stream2;
  file_stream2.open(oupFileDisp.c_str(),std::ios::out);
  ////2.2.1) dump the file infomation
  
  file_stream2 << d_num_axial << "\t" << 3 
    << "  #" << "average disp, min disp, max disp"
    << "\n";
  //// 2.2.2) dump all the data on axial displacement
    for (it_level=0; it_level<total_levels; it_level++)
    {
      
      int it_para;
      for (it_para=0; it_para< 3; it_para++)
      {file_stream2<< vec3_axial_length[it_level][it_para] << "\t";
      }
      file_stream2 << "\n";
	
    }
  file_stream2.close();
  std::cout << OutPrefix<<"Done: write axial disp. in file:" << oupFileDisp<< std::endl;
  ////2.1.2) axial spring_length << 05/2/4/2013
  ////2.2) start write the value
  std::cout << OutPrefix<<"Begin: write area in file:" << oupFileArea << std::endl;
  std::ofstream file_stream;
  file_stream.open(oupFileArea.c_str(),std::ios::out);
  ////2.2.1) dump the file infomation
  
  file_stream << d_num_axial << "\t" << num_para+1 
    << "  #" << d_fcn_info_map[d_fcn_id]
    << "\n";
  //// 2.2.2) dump all the area and parameters
    for (it_level=0; it_level<total_levels; it_level++)
    {
      
      int it_para;
      for (it_para=0; it_para< num_para; it_para++)
      {file_stream<< vector_para[it_level][it_para] << "\t";
      }
      file_stream << "\n";
	
    }
 //// 2.2.3) finish all the dump
  file_stream.close();
  std::cout << OutPrefix<<"Done: write area in file:" << oupFileArea << std::endl;
  //vector_para.clear();
  vector_pts.clear();
  // 2.2.4) calculate the volume
  std::cout<< OutPrefix<<"Begin: calculate volume:fcn_id:" << d_fcn_id << std::endl;
  //volume=0.0;
  (d_volume_fcn_map[d_fcn_id])(vector_para,total_levels,volume[ind_layer]);
  std::cout << OutPrefix<<"Done: calculate volume" << std::endl;
  // record the area and axial rest_length;
  for(it_level=0; it_level<total_levels; it_level++)
  {vec_summary[it_level][ind_layer*2]=vector_para[it_level][0]; // area
   vec_summary[it_level][ind_layer*2+1]=vec3_axial_length[it_level][0];   // axial displacement
  }
  vec3_axial_length.clear();
  vector_para.clear();
  } //finish loop ind_layer
  //// 2.4) dump the summary file;
  std::cout << OutPrefix<<"Begin: write summary in file:" << oupFileSum << std::endl;
  std::ofstream file_stream4;
  file_stream4.open(oupFileSum.c_str(),std::ios::out);
  ////2.2.1) dump the file infomation
  int num_para=2*d_num_layers;
  file_stream4 << d_num_axial << "\t" << d_num_layers 
    << "  #" << "num_levels, num_layers (1.area,aixal 2.displacement"
    << "\n";
  //// 2.2.2) dump all the area and displacement
    int it_level;
    for (it_level=0; it_level<total_levels; it_level++)
    {
      
      int it_para;
      for (it_para=0; it_para< num_para; it_para++)
      {file_stream4<< vec_summary[it_level][it_para] << "\t";
      }
      file_stream4 << "\n";
	
    }
 //// 2.2.3) finish all the dump
  file_stream4.close();
  std::cout << OutPrefix<<"Done: write summary in file:" << oupFileSum << std::endl;
  vec_summary.clear();
  IsWriteDone=true;
 // std::cout << OutPrefix<<"Done: IsWriteDone= "<< IsWriteDone << std::endl;
  return IsWriteDone;

}// Caculate and dump area function ptr to calculate the area




