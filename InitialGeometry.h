// generate xda file for the mesh
// created by 01/23/2014 walter
// input: infor: nodal, element, boundary conditions
// output: generate filename.xda
#ifndef _INCL_INIT_GEOMETRY
#define _INCL_INIT_GEOMETRY
// C++ includes
#include <iostream>
#include <string>
#include "CPts.h"
#include <vector>

using namespace std;
namespace LibGeometry {


static const double PI = 3.1415926;
static const string ErrorMsg= "[Error:]";
static const string WarnMsg= "[Warning:]";


// **** structure of connectivity
struct ConnectInfo
{ vector <int> node_ids;
  int no_nodes; //number of nodes
  int elem_id;
  int pare_elem_id; // default -1

};
// multiple definition of functions in *.h file---> use inline to avoid
// read from input stream
inline
 istream& operator>> (istream& is, ConnectInfo & connect_info)
 {  cout << "not implemented" << endl;

    return is;
}
// multiple definition of functions in *.h file---> declare in h file, define it in cpp file
// write to output stream
inline
ostream& operator<<(ostream& os, const ConnectInfo & connect_info)
{ if (connect_info.no_nodes < 1) cout << ErrorMsg << "empty connectitvity info!" << endl;
    else{ // output node id
        for (int kid =0; kid < connect_info.no_nodes; kid ++)
        { os << connect_info.node_ids[kid] << "\t";
        }
        os << connect_info.elem_id << "\t";
        os << connect_info.pare_elem_id << "\t"; // still on that line;

    }

  return os;

}


// **** structure of blocked connectivity
struct BlockInfo
{ int no_patches; // no. of the element in this block
  vector<ConnectInfo> connect_info_data; // number of connectivity info
  int level_id; // = 0; // level id, defalt the 0 level
  string block_name ;
  int patch_type_id; // patch type id;
};

// read from input stream
inline
 istream& operator>> (istream& is, BlockInfo & block_info)
 {  cout << "not implemented" << endl;

    return is;
}
// write to output stream
inline
 ostream& operator<<(ostream& os,  const BlockInfo & block_info)
{ if (block_info.no_patches<1) cout << ErrorMsg << "empty block info!" << endl;
  else{ // output connect info
        for (int kid =0; kid < block_info.no_patches; kid++)
        {   os << (block_info.connect_info_data)[kid];
            os << "#:" << block_info.patch_type_id << "\t" << block_info.block_name << "\n";
        }


    }

  return os;
}
// *** structure of boundary condition information
// 3 short integers: elem_id, face_id, boundary_id
 struct BCInfo {
     int elem_type_id ;//=0;
     int face_id ;//=0;
     int bc_id ;//=0;
 };
 inline
 ostream& operator<<(ostream& os, const BCInfo & bc_info)
 { os << bc_info.elem_type_id << "\t" << bc_info.face_id << "\t" << bc_info.bc_id << "\t";

 }
// *** structure of layer mesh information
// consistent with previous design: layer by layer
// notice here, layer = block; and we use layer to compute mesh info on each layer
// default 3 dimension problem
struct LayerInfo {
    int cs_type_id ; //=0; // coordinate system type, 0: cylindrical type
    int grid_info_method ; //=0; // info on grid: 0: gridNum, meshSize;
    double vec_size[6]; // size of the layer; r, z, theta;
    double vec_meshSize[3]; // ds1, ds2, ds3, we get the maixmal value, if vec_grid NUm is given
    int vec_gridNum[3]; // Ns1, Ns2, Ns3
    int patch_type_id; // patch-type, by default = element type in block info
    int patch_size; // number of nodes on each patch, by default = number of nodes for each element type
};

class GenerateMesh
{
public:
    // constructor
    GenerateMesh(string InputFile);
    // read inputInfo
    bool  ReadInfoFromFile(string InputFile);
    // output info
    bool GenerateXDAFile(); // this is the whole XDA file for the whole structure
    // check info

    string XDAFileName(){return d_op_filename + ".xda";}

    // print overall info
    void print_info();
    void print_elem_info();
    void print_node_info();
    int get_num_parts()
    {return d_no_block;}
    int get_no_dim()
    {return d_n_dim;}

private:
    // private method
    bool build_nodes_per_part(vector<MY_APP::CPts>& vec_local_nodes, int part_id, int no_nodes_skipped = 0);
    // build elements (patches or nodal connectivity) based on layer_info
    // we first adopt local node index/element index: node start =0; element_start =0;
    bool build_elements_per_part(vector<ConnectInfo>& vec_local_patches, int part_id, int no_nodes_offset =0, int no_elements_offset =0);
    // overall info
    int d_n_dim ; //=3;
    string d_op_filename ;//="mesh";


    // overall info on infoFile;
    int d_version_no ;//=1;
    int d_noUsed_lines; // =0;
    // header info- part0;
    int d_refine_levels; // = 0; // by default, we do not use refinement
    int d_no_elem;
    int d_no_nodes;

    int d_no_bcs;


    // data info- part 1: connectivity
    // each line of connectivity info
    // node ids; elem_id; paranet_elem_id

    vector<LayerInfo> d_layers;
    int d_no_layer; // always number of layers = d_no_block

    vector<BlockInfo> d_blocks;
    int d_no_block; // number of blocks
    int d_no_level; //=1; // number of levels

    // data info- part 2: node coordinates
    int d_num_node;
    vector<MY_APP::CPts> d_vec_node;

    // data info- part 3: BC info


    int d_num_bcs;
    vector<BCInfo> d_vec_bc;
    // constructor
    GenerateMesh();
    GenerateMesh(GenerateMesh&);
    //~GenerateMesh();
    GenerateMesh & operator = (GenerateMesh &);

};// class GenerateMesh

}
#endif
