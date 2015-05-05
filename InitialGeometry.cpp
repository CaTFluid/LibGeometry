// generate xda file for the mesh
// created by 01/23/2014 walter
// input: infor: nodal, element, boundary conditions
// output: generate filename.xda
#include "DiscardComments.h"
//#include "DiscardComments.h"
#include "InitialGeometry.h"
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
using namespace std;
using namespace LibGeometry;
using namespace MY_APP;
using namespace MY_CPP;



// for classes: GenerateMesh
bool  GenerateMesh::ReadInfoFromFile(string InfoFile)
{

    std::string line_string;
    std::ifstream file_stream;
    std::string FileName=InfoFile;
    std::string HeaderPrefix = "[HEARER_INFO]:";
    std::string PartPrefix="[BLOCK_INFO]:";
    std::string InfoPrefix="[NOUSED_INFO]:";

//    std::string STRING_CHECK_READ= "\t+++[Check_READ]:";
//    std::string STRING_CHECK_IMPORTANT="\t+++";
//    std::string STRING_DIVIDE_LINE(50,'*');

    file_stream.open(FileName.c_str(),std::ios::in);//input(to the screen)=read
    if (file_stream.is_open())
   {  std::istringstream line_stream;
        bool dummy_bool =true;

        // 0) first lines: version number + info lines
        if(!std::getline(file_stream,line_string))
            // error
        {std::cout<< STRING_ERROR<<" for read"<<std::endl;

            return false;
        }
        else
        {
            std::cout<< InfoPrefix << line_string << std::endl;
            line_string=discard_comments(line_string);
            line_stream.str(line_string);
            line_stream >> this->d_version_no;
            line_stream >> this->d_noUsed_lines;
            for (int k_line =0; k_line < this->d_noUsed_lines; k_line ++)
            {   if (std::getline(file_stream,line_string))
                    std::cout<< InfoPrefix << line_string<< std::endl;
            }

        }
        // 1) header info: info on output file, dimension
        if(!std::getline(file_stream,line_string))
            // error
        {std::cout<< STRING_ERROR <<  HeaderPrefix<<" on output file, dimension number-- error for read" <<std::endl;
            return false;
        }
        else
        {   //  title info
            std::cout<< HeaderPrefix << line_string << std::endl;
            // output file
            dummy_bool = std::getline(file_stream,line_string);
            std::cout<< HeaderPrefix << line_string << std::endl;
            line_string=discard_comments(line_string);
            line_stream.str(line_string);
            line_stream >> this->d_op_filename;
            std::cout<< STRING_CHECK_READ << "output XDA file:" << this->d_op_filename << ".xda" << std::endl;
            // dimension
            dummy_bool = std::getline(file_stream,line_string);
            std::cout<< HeaderPrefix << line_string << std::endl;
            line_string=discard_comments(line_string);
            line_stream.str(line_string);
            line_stream >> this->d_n_dim;
            std::cout<< STRING_CHECK_READ << "dimension:" << this->d_n_dim << std::endl;


        }

        std::cout << STRING_DIVIDE_LINE << endl;
        // 2) data part: info on Block/Layers
        if(!std::getline(file_stream,line_string))
            // error
        {std::cout<< PartPrefix <<" on block/layer-- error for read" <<std::endl;
            return false;
        }
        else
        {   // part title info
            std::cout<< PartPrefix << line_string << std::endl;
            // data
            dummy_bool = std::getline(file_stream,line_string);
            line_string=discard_comments(line_string);
            line_stream.str(line_string);
            line_stream >> this->d_no_block;
            this->d_blocks.resize(this->d_no_block);
            this->d_no_layer=this->d_no_block;
            this->d_layers.resize(this->d_no_layer);
            std::cout<<PartPrefix  << "d_no_layer =" << this->d_no_layer << std::endl;
           line_stream.clear();
            // read info for each layer
            for (int it_layer =0 ; it_layer<this->d_no_layer; it_layer ++)
            {
                std::cout << STRING_DIVIDE_LINE << endl;
                // print layer id info
                LayerInfo* a_layer_info = &(this->d_layers[it_layer]);
                dummy_bool = std::getline(file_stream,line_string);
                std::cout<<"+++++ Layer id= "<< it_layer << endl;
                std::cout<<"+++++ check header:" << line_string << endl;
                // put it into the structure
                dummy_bool = std::getline(file_stream,line_string);
                std::cout<< PartPrefix << line_string << std::endl;
                line_string=discard_comments(line_string);
                line_stream.str(line_string);
                line_stream >> a_layer_info->cs_type_id
                            >> a_layer_info->grid_info_method;
                line_stream.clear();
                if (a_layer_info->grid_info_method > 0)
                    std::cout<< PartPrefix  << "not implemented except given grid info" << std::endl;
                else{
                    // domain size
                    dummy_bool = std::getline(file_stream,line_string);
                    std::cout<< PartPrefix << line_string << std::endl;
                    line_string=discard_comments(line_string);
                    line_stream.str(line_string);
                    line_stream >> a_layer_info->vec_size[0]  // r_i
                                >> a_layer_info->vec_size[1]  // r_o
                                >> a_layer_info->vec_size[2]  // z_lower
                                >> a_layer_info->vec_size[3]  // z_upper
                                >> a_layer_info->vec_size[4]  // not used theta_0 =0
                                >> a_layer_info->vec_size[5]; // not used theta-2pi = 2pi
                    std::cout<< STRING_CHECK_READ;
                    line_stream.clear();
                    a_layer_info->vec_size[5] = 2*PI;
                    for (int it_size =0 ; it_size < 6; it_size ++)
                        std::cout << "s" << it_size << "=" << a_layer_info->vec_size[it_size]
                                  << "; ";
                    std::cout<<std::endl;
                    // grid size
                    dummy_bool = std::getline(file_stream,line_string);
                    std::cout<< PartPrefix << line_string << std::endl;
                    line_string=discard_comments(line_string);
                    std::cout<< PartPrefix << " discard :"<< line_string << std::endl;
                    line_stream.str(line_string);
                    line_stream >> a_layer_info->vec_gridNum[0]
                                >> a_layer_info->vec_gridNum[1]
                                >> a_layer_info->vec_gridNum[2];

                    std::cout<< STRING_CHECK_READ;
                    line_stream.clear();
                    for (int it_size =0 ; it_size < 3; it_size ++)
                        std::cout << "Ns" << it_size << "=" << a_layer_info->vec_gridNum[it_size]
                                  << "; ";
                    std::cout<< std::endl;
                    // dr:
                    a_layer_info->vec_meshSize[0]=(a_layer_info->vec_size[1]
                                                   - a_layer_info->vec_size[0])/
                            ( a_layer_info->vec_gridNum[0]);
                    // dz:
                    a_layer_info->vec_meshSize[1]=(a_layer_info->vec_size[3]
                                                   - a_layer_info->vec_size[2])/
                            ( a_layer_info->vec_gridNum[1]);
                    // rdtheta:
                    a_layer_info->vec_meshSize[2]=(a_layer_info->vec_size[1])* 2.0 * 3.1416 /
                            ( a_layer_info->vec_gridNum[2]);
                    std::cout << PartPrefix<< "mesh size (computed)" << std::endl;
                    std::cout<< STRING_CHECK_READ;
                    for (int it_size =0 ; it_size < 3; it_size ++)
                        std::cout << "Ds" << it_size << "=" << a_layer_info->vec_meshSize[it_size]
                                  << "; ";
                    std::cout << std::endl;
                    // element/patch info: patch type id; connected nodes number (optional)
                    dummy_bool = std::getline(file_stream,line_string);
                    std::cout<< PartPrefix << line_string << std::endl;
                    line_string=discard_comments(line_string);
                    line_stream.str(line_string);
                    line_stream >> a_layer_info->patch_type_id
                                >> a_layer_info->patch_size;
                    line_stream.clear();

                }

            }



        }

        // 3) finish file
        file_stream.close();
        cout<<PartPrefix  << "  success" << endl;
        return true;

    }
    return false;
} // GenerateMesh::ReadInfoFromFile(string InfoFile)
// with read info [d_layers], we generate nodes part per part
bool  GenerateMesh::build_nodes_per_part(std::vector<CPts>& vec_local_nodes, int part_id, int no_nodes_skipped)
{// here each part = each layer

    LayerInfo * a_layer_info = &(this->d_layers[part_id]);
    if (a_layer_info->cs_type_id !=0 )
    {    std::cout << "[Error]:not cylindrical cs" << std::endl;
        return false;
    }
    else
    {
        std::cout << STRING_CHECK_IMPORTANT << "Only support uniform dr, dz, dtheta" << std::endl;
        int Nr = a_layer_info->vec_gridNum[0];
        int Nz = a_layer_info->vec_gridNum[1];
        int Ntheta = a_layer_info->vec_gridNum[2];

        double ri=a_layer_info->vec_size[0];
        double ro=a_layer_info->vec_size[1];
        double z1=a_layer_info->vec_size[2];
        double z2=a_layer_info->vec_size[3];
        double t1=0.0;
        double t2= PI * 2;
        //
        int total_nodes = (Nr+1) * (Nz+1) * Ntheta;
        double dr= (ro-ri)/Nr;
        double dz= (z2-z1)/Nz;
        double dtheta = (t2 -t1) / Ntheta;
        double r_start = ri;
        int Nr_start =0;
        if (no_nodes_skipped >0)
        {

            if (Ntheta* (Nz+1) - no_nodes_skipped !=0)
            {
            std::cout << STRING_WARN << "Only support: skip the first layer of nodes" << std::endl;
            return false;
            }
            r_start = ri; // + dr; // r_cod for the firs new node; since it_layer =1 account for this
             Nr_start =1; // start from the second layer;

        }
        // build nodes
        vec_local_nodes.resize(total_nodes - no_nodes_skipped);
        // follow before logic, we loop r, loop theta, loop z;
        // loop r_cods first;
        int node_index =-1;
        for (int it_layer = Nr_start; it_layer<(Nr+1); it_layer++)
        {
            double r_cod=r_start + (it_layer*dr);
            // loop t_cod second
            for (int it_theta=0; it_theta<Ntheta; it_theta++)
            {
                double t_cod = it_theta * dtheta;
                double x_cod = cos(t_cod) * r_cod;  //sin(t_cod) * r_cod;
                double y_cod= sin(t_cod) * r_cod;//cos(t_cod) * r_cod;
                // loop z_cod third;
                for (int it_level = 0; it_level<(Nz+1); it_level++)

                {   node_index++;
                    double z_cod = it_level * dz + z1;
                    vec_local_nodes[node_index].x=x_cod;
                    vec_local_nodes[node_index].y=y_cod;
                    vec_local_nodes[node_index].z=z_cod;
                }
            }
        }// for layer
    return true;

    }

} // GenerateMesh::build_nodes_per_part()

bool GenerateMesh::build_elements_per_part(vector<ConnectInfo>& vec_local_patches, int part_id, int no_nodes_offset, int no_elem_offset )

{ // the label of elments on some part -> subdomian Id of the elements
  // part 1) analyze the layer info

    LayerInfo * a_layer_info = &(this->d_layers[part_id]);
    if (a_layer_info->cs_type_id !=0 )
    {    std::cout << "[Error]: not cylindrical cs" << std::endl;
        return false;
    }
    std::cout <<STRING_CHECK_IMPORTANT<< "Only support uniform dr, dz, dtheta" << std::endl;
    int Nr = a_layer_info->vec_gridNum[0];
    int Nz = a_layer_info->vec_gridNum[1];
    int Ntheta = a_layer_info->vec_gridNum[2];

//    double ri=a_layer_info->vec_size[0];
//    double ro=a_layer_info->vec_size[1];
//    double z1=a_layer_info->vec_size[2];
//    double z2=a_layer_info->vec_size[3];
//    double t1=0.0;
//    double t2= PI * 2;
    int num_nodes = a_layer_info->patch_size;
    if (num_nodes <1)
    { std::cout <<STRING_WARN<< "use default hex8 for 3D problem" << std::endl;
        num_nodes = 8;
    }
    if ( this->d_n_dim ==3) //
    {
        int total_elements = Nr * Nz * Ntheta;
        vec_local_patches.resize(total_elements);
        // >> version 2 modified 07/24/2014
	//int elem_index =-1;
	int elem_index =-1;
	// << version 2 modified 07/24/2014
        int node_start =no_nodes_offset;
        // loop r
        for (int it_layer =0 ; it_layer < Nr; it_layer ++)
        { int node_layer = (it_layer * (Nz+1) * Ntheta );
            for (int it_theta =0; it_theta<Ntheta; it_theta ++)
            {  int  node_theta= (it_theta *(Nz+1));
                // notice, Ntheta -> 0 connected theta = 0;
                int next_theta = ((it_theta+1)==Ntheta? 0: it_theta+1);
                // node (r, z, t2) - node (r, z, t1)
                int id_diff= (next_theta - it_theta) * (Nz+1);

                for (int it_level =0; it_level<Nz; it_level ++)
                { elem_index++;
                    ConnectInfo* a_element = &(vec_local_patches[elem_index]);
                  a_element->no_nodes=num_nodes;
                  a_element->node_ids.resize(num_nodes);
                  a_element->elem_id = elem_index + no_elem_offset; // local elem_id
                  a_element->pare_elem_id = -1; //( default no parent element)
                  int node_first = node_start + node_layer+node_theta+(it_level); // inner_lower_left node index;
                  // for hex8
                  a_element->node_ids[0] = node_first; // (r1, z1, t1)
                  a_element->node_ids[1] = node_first + (Nz+1) * Ntheta; // (r2, z1, t1)
                  a_element->node_ids[2] = node_first + (Nz+1) * Ntheta +id_diff; // (r2, z1,t2)
                  a_element->node_ids[3] = node_first + id_diff; // (r1, z1, t2)
                  a_element->node_ids[4] = node_first + 1; // (r1, z2, t1)
                  a_element->node_ids[5] = node_first + (Nz+1) * Ntheta +1 ;  // (r2, z2, t1)
                  a_element->node_ids[6] = node_first + (Nz+1) * Ntheta +id_diff + 1;  // (r2, z2, t2)
                  a_element->node_ids[7] = node_first + id_diff +1; // (r1,z2,t2)
                }
            }
        }


    }
    else // other type id
    {
      std::cout << STRING_WARN<<"Only support 3 dim Hex8" << std::endl;
      return false;
    }
    
    return true;

} // GenerateMesh::build_elements_per_part()

// constructor to obtain the data from inputFile
GenerateMesh::GenerateMesh(string InputFile)
    :d_n_dim(3)
    ,d_no_level(1)
    ,d_version_no(1)
    ,d_noUsed_lines (0)
    ,d_refine_levels (0)
    ,d_op_filename ("mesh")
{

    if (!GenerateMesh::ReadInfoFromFile(InputFile))
    {
        std::cout <<STRING_ERROR << " fail to read file" << std::endl;
        return;
    }

    std::cout <<"[Begin]:   build/write separate nodes/elements data" << endl;
    for (int part_id =0 ; part_id<this->d_no_layer; part_id++)
    {   // output file info
        char strlayer[10];
        sprintf(strlayer,"%01d", part_id);
        string part_fix= "_part";
        string part_file_name = this -> d_op_filename + part_fix + strlayer;
        // node info
        vector<CPts> vec_local_nodes;
        if (this->build_nodes_per_part(vec_local_nodes,part_id,0))
        { // dump to the file:
            std::cout << STRING_WRITE_BEGIN<<" write nodes part " << part_id << " of " << d_no_layer << endl;
            int node_size = vec_local_nodes.size();
            string node_file = part_file_name + ".node";
            std::ofstream file_stream2;
            file_stream2.open(node_file.c_str(),std::ios::out);
            // header info
            file_stream2 << node_size << "\t" << 3
              << "  #" << "x, y , z"
              << "\n";
            for (int it_node =0 ; it_node < node_size ; it_node++)
            {
                file_stream2 << vec_local_nodes[it_node] << "\n";
            }
            file_stream2.close();
            std::cout << STRING_WRITE_END<< "write nodes part " << part_id << " of " << d_no_layer << endl;
        }
        // element info
        vector<ConnectInfo> vec_local_patches;

        if (this->build_elements_per_part(vec_local_patches, part_id, 0))
        { // dump to the file:
            int patch_size = vec_local_patches.size();
            //cout << "Check patch size: " << patch_size << endl;
            std::cout <<  STRING_WRITE_BEGIN<<" write elements part " << part_id << " of " << d_no_layer << endl;

            string patch_file = part_file_name + ".element";
            std::ofstream file_stream2;
            file_stream2.open(patch_file.c_str(),std::ios::out);
            // header info
            file_stream2 << patch_size << "\t" << vec_local_patches[0].no_nodes + 2
              << "  #" << "connected node index, elem_id, pare_elem_id = -1"
              << "\n";
            for (int it_patch=0; it_patch < patch_size; it_patch++)
            {
                file_stream2 << vec_local_patches[it_patch] << "\n";
            }
            file_stream2.close();
            std::cout <<STRING_WRITE_END <<" write elements part " << part_id << " of " << d_no_layer << endl;
        }
        else
        {
            std::cout << STRING_ERROR<<"fail to build element?" << endl;
        }
        // xda file info:
        // currently, we ignore the BC info
        if (vec_local_nodes.size()>0 && vec_local_patches.size()>0)
        {
            std::cout << STRING_WRITE_BEGIN <<"write elements part " << part_id << " of " << d_no_layer << endl;
            string xda_file = part_file_name+".xda";
            std::ofstream file_stream2;
            file_stream2.open(xda_file.c_str(),std::ios::out);
            // header info
            file_stream2 << "LIBM" << "\t" << "0" << "\n";
            file_stream2 << vec_local_patches.size() << "\t" << "# Num. Elements" << "\n";
            file_stream2 << vec_local_nodes.size() << "\t" << "# Num. Nodes" << "\n";
            file_stream2 << (vec_local_patches[0].no_nodes +2 )* vec_local_patches.size() << "\t" << "# Length of connectivity vector" << "\n";
            file_stream2 << "0" << "\t" << "# Num. Boundary Conds." << "\n";
            file_stream2 << "65536\t# String size (igonore)" << "\n";
            file_stream2 << "1" << "\t # Num. Element Blocks"<< "\n";
            file_stream2 << this->d_layers[part_id].patch_type_id << "\t # Element types in each block"<< "\n";
            file_stream2 << vec_local_patches.size() << "\t # Num. of elements in each block at each level"<< "\n";
            file_stream2 << "Id String" << "\n";
            file_stream2 << "Title String" << "\n";
            // data info: connectivity
            int patch_size = vec_local_patches.size();
            int node_size = vec_local_nodes.size();
            for (int it_patch=0; it_patch < patch_size; it_patch++)
            {
                file_stream2 << vec_local_patches[it_patch] << "\n";
            }
            // data info: nodes
            for (int it_node =0 ; it_node < node_size ; it_node++)
            {
                file_stream2 << vec_local_nodes[it_node] << "\n";
            }
            file_stream2.close();

        }
        else
        {
            std::cout << STRING_ERROR << " fail to obtain node or element info" << std::endl;
        }


    }
    std::cout <<STRING_WRITE_END <<" build/write separate nodes/elements data" << endl;
    std::cout << STRING_DIVIDE_LINE << endl;
    //std::cout <<"[Begin]:   build/write the whole structure" << endl;
}// constructor

// generate the XDA file for the whole structure: 
 bool GenerateMesh::GenerateXDAFile()
 {// input_info
 std::cout << STRING_DIVIDE_LINE << endl;
 std::cout <<"[Begin]:   build/write collected nodes/elements data" << endl;
 int num_parts = this->d_no_layer;
 std::vector<int> vec_layer_offset; // for the element; change node_id;
 std::vector<int> vec_layer_skipped_nodes; // for the nodes
 std::vector<int> vec_layer_offset_element; // for the element; change element_id
 
 std::vector< std::vector<CPts> > block_vec_nodes;
 std::vector< std::vector<ConnectInfo> > block_vec_patches;
 
 // deal with node id consistency
 vec_layer_skipped_nodes.resize(num_parts); 
 vec_layer_offset.resize(num_parts);
 vec_layer_offset_element.resize(num_parts);
 
 block_vec_nodes.resize(num_parts);
 block_vec_patches.resize(num_parts);
 
 
 
 unsigned int sum_node =0;
 unsigned int sum_elem=0;
 unsigned int curPart_node_id_begin =0;
 unsigned int curPart_node_id_end =0;
 unsigned int curPart_element_id_begin =0;
 for (int part_id =0 ; part_id<num_parts; part_id++)
 {
   // current part node info;
   LayerInfo * cur_layer_info = &(this->d_layers[part_id]); // last layer info
             int Nr = cur_layer_info->vec_gridNum[0];
        int Nz = cur_layer_info->vec_gridNum[1];
        int Ntheta = cur_layer_info->vec_gridNum[2];
	int CurLayerTotalNodes = (Nr +1 ) * (Nz +1) * Ntheta;
	std::cout << STRING_DIVIDE_LINE << endl;
   
   
   if (part_id ==0 ) // first layer
   { // no offset of node id
    curPart_node_id_begin =0; // this is for patch info; we consider all the layers
    curPart_node_id_end =curPart_node_id_begin + CurLayerTotalNodes -1;
    curPart_element_id_begin =0;
    
    vec_layer_skipped_nodes[part_id]=0; // do not skippped nodes    
    vec_layer_offset[part_id]=curPart_node_id_begin;
    vec_layer_offset_element[part_id]=curPart_element_id_begin;
   }
   else
   { LayerInfo * a_layer_info = &(this->d_layers[part_id -1]); // last layer info
             int a_Nr = a_layer_info->vec_gridNum[0];
        int a_Nz = a_layer_info->vec_gridNum[1];
        int a_Ntheta = a_layer_info->vec_gridNum[2];
	int LastLayerTotalNodes = (a_Nr +1 ) * (a_Nz +1) * a_Ntheta;
	int LastSurfaceNodes = (a_Nz+1) * a_Ntheta;
	// info from last layer
	curPart_node_id_begin = curPart_node_id_end - LastSurfaceNodes+1;
	curPart_node_id_end =curPart_node_id_begin + CurLayerTotalNodes -1;
	curPart_element_id_begin =curPart_element_id_begin + block_vec_patches[part_id-1].size();
	
	vec_layer_offset[part_id] = curPart_node_id_begin; // for patch info: this node 0 --> changed
	vec_layer_offset_element[part_id] = curPart_element_id_begin; // for patch info: element 0 --> changed
	vec_layer_skipped_nodes[part_id]=LastSurfaceNodes; // for collecting nodes: skip repeated nodes
 
   }
   
   
   // begin to collect all the nodes
   if (!this->build_nodes_per_part(block_vec_nodes[part_id],part_id,vec_layer_skipped_nodes[part_id]))
   { std::cout << STRING_ERROR<<"fail to build nodes from part #" << part_id << endl;
   }
   else
   { std::cout << STRING_CHECK_IMPORTANT << "From part #" << part_id <<", get: "<< block_vec_nodes[part_id].size() <<" nodes, and skip " << vec_layer_skipped_nodes[part_id] << " nodes" << endl; 
   sum_node = sum_node + block_vec_nodes[part_id].size();
   }
   
   if (!this->build_elements_per_part(block_vec_patches[part_id], part_id, vec_layer_offset[part_id],vec_layer_offset_element[part_id]))
   { std::cout << STRING_ERROR<<"fail to build elements from part #" << part_id << endl;
   }
   else
   {std::cout << STRING_CHECK_IMPORTANT << "From part #" << part_id <<", get: "<< block_vec_patches[part_id].size() <<" elements, and node id starts with " << vec_layer_offset[part_id] <<  endl; 
   sum_elem =sum_elem +block_vec_patches[part_id].size(); 
   
   }
    std::cout << STRING_CHECK_IMPORTANT <<  "part #: "<<part_id <<"; node_id_start, nodes_total, node_id_end: " 
	<< curPart_node_id_begin <<"; "<<CurLayerTotalNodes << "; " 
	<< curPart_node_id_end << std::endl;
   
	
    std::cout << STRING_CHECK_IMPORTANT <<  "part #: "<<part_id <<"; elem_id_start,  elem_id_end: " 
	<< curPart_element_id_begin <<"; "<<curPart_element_id_begin + block_vec_patches[part_id].size() -1 << std::endl;	
	
   
 
   
 } // loop part_id;
 
 // 2) begin to write all the data:
 std::cout << STRING_WRITE_BEGIN<< "write all the nodes " << endl;
 string node_file =this -> d_op_filename +".node";
 std::ofstream file_stream2;
            file_stream2.open(node_file.c_str(),std::ios::out);
            // header info
            file_stream2 << sum_node << "\t" << 3
              << "  #" << "x, y , z"
              << "\n";
	      
for (int part_id =0 ; part_id<num_parts; part_id++)
{
  for (int it_node =0 ; it_node < block_vec_nodes[part_id].size() ; it_node++)
            {
                file_stream2 << block_vec_nodes[part_id][it_node] << "\n";
            }
}

file_stream2.close();
std::cout << STRING_WRITE_END<< "write all the nodes " << endl;
 
 // 3) begin to write all the elements
	    
 // begin to write all the data:
 std::cout << STRING_WRITE_BEGIN<< "write all the elements " << endl;
 string patch_file =this -> d_op_filename +".element";
 //std::ofstream file_stream2;
            file_stream2.open(patch_file.c_str(),std::ios::out);
            // header info
            file_stream2 << sum_elem<< "\t" << block_vec_patches[0][0].no_nodes + 2
              << "  #" << "connected node index, elem_id, pare_elem_id = -1"
              << "\n"; 
for (int part_id =0 ; part_id<num_parts; part_id++)
{
  for (int it_patch=0; it_patch < block_vec_patches[part_id].size(); it_patch++)
  {
    file_stream2 << block_vec_patches[part_id][it_patch] << "\n";
  }
}

file_stream2.close();
std::cout <<STRING_WRITE_END <<" write all the elements " << endl;
 
// 4) begin to write xda_file
std::cout << STRING_WRITE_BEGIN<< "write XDA info " << endl;

string xda_file = this -> d_op_filename+".xda";
            
            file_stream2.open(xda_file.c_str(),std::ios::out);
            // header info
            file_stream2 << "LIBM" << "\t" << "0" << "\n";
            file_stream2 << sum_elem << "\t" << "# Num. Elements" << "\n";
            file_stream2 << sum_node << "\t" << "# Num. Nodes" << "\n";
            file_stream2 << ( block_vec_patches[0][0].no_nodes+2 )* sum_elem << "\t" << "# Length of connectivity vector" << "\n";
            file_stream2 << "0" << "\t" << "# Num. Boundary Conds." << "\n";
            file_stream2 << "65536\t# String size (igonore)" << "\n";
            file_stream2 << "1" << "\t # Num. Element Blocks"<< "\n";
            file_stream2 << this->d_layers[0].patch_type_id << "\t # Element types in each block"<< "\n";
            file_stream2 << sum_elem << "\t # Num. of elements in each block at each level"<< "\n";
            file_stream2 << "Id String" << "\n";
            file_stream2 << "Title String" << "\n";
            // data info: connectivity
for (int part_id =0 ; part_id<num_parts; part_id++)
{
  for (int it_patch=0; it_patch < block_vec_patches[part_id].size(); it_patch++)
  {
    file_stream2 << block_vec_patches[part_id][it_patch] << "\n";
  }
}
// data info: nodes

for (int part_id =0 ; part_id<num_parts; part_id++)
{
  for (int it_node =0 ; it_node < block_vec_nodes[part_id].size() ; it_node++)
            {
                file_stream2 << block_vec_nodes[part_id][it_node] << "\n";
            }
}
file_stream2.close();
std::cout <<STRING_WRITE_END <<" write XDA info " << endl;
std::cout <<STRING_WARN <<" +++++++++++ for material.template +++++++++++++ " << endl;
std::cout <<STRING_WARN <<" +++++++++++ Note: info. for material.template +++++++++++++ " << endl;

 string log_file =this -> d_op_filename +".log";
 std::ofstream file_stream_log;
            file_stream_log.open(log_file.c_str(),std::ios::out);

std::cout <<STRING_WARN <<" +++++++++++ is also written into the log file  +++++++++++++ " << endl;

// >> first print out the overall info. of the geometry: needed for geometry.info
file_stream_log<< "# version 3" << ": 3D multiple layered cylindrical tube " << "\n";
file_stream_log<< "please modify geometry(material).info, manually or with script" << "\n";
// 1) part 1: begin boundary_info;
double vec_tube_size[4];
int last_layer = this->d_no_layer - 1;
vec_tube_size[0] = this->d_layers[0].vec_size[0]; // r_i
vec_tube_size[1] = this->d_layers[last_layer].vec_size[1]; //r_o
vec_tube_size[2] = this->d_layers[0].vec_size[2]; // z_lower
vec_tube_size[3] =this->d_layers[0].vec_size[3]; // z_upper
file_stream_log<< "#// top, bot, r_inner, r_outer  3D multiple layered cylindrical tube" << "\n";
file_stream_log<< vec_tube_size[3] << "\t" << vec_tube_size[2] 
	       << "\t" << vec_tube_size[0] << "\t" << vec_tube_size[1] << "\n";
// 2) part 2: star_id;
// 
std::vector<int> vec_start_ids(num_parts);
vec_start_ids[0]= 0;
file_stream_log << num_parts << "\t #number of domains" << "\n";


int pre_elem_sum =0;
cout << "block/layer # 0" << ": elem_id_bgn, elem_id_end: " << pre_elem_sum << "\t" <<
block_vec_patches[0].size() + pre_elem_sum -1<< endl;
file_stream_log <<  vec_start_ids[0] << "\t";


// for material infor
for (int part_id =1 ; part_id<num_parts; part_id++)
{
  pre_elem_sum = pre_elem_sum + block_vec_patches[part_id-1].size();
cout << "block/layer #" << part_id << ": elem_id_bgn, elem_id_end: " << pre_elem_sum << "\t" <<
block_vec_patches[part_id].size() + pre_elem_sum -1<< endl;
// >> version 3 also cout to log file
vec_start_ids[part_id]= pre_elem_sum;
file_stream_log <<  vec_start_ids[part_id] << "\t";

}
file_stream_log<< "#start_id;" << "\n";
// 3) we add the last element id;
file_stream_log << pre_elem_sum + block_vec_patches[num_parts-1].size() -1 << "\t #end_id of the whole structure;" << "\n"; 


// 4) for plot to check the geometry (needs python)
file_stream_log <<"layer_id, r_i, r_o, nr, ntheta" 
		<<";with L and dz: " << this->d_layers[0].vec_size[3] 
		<<" and " << this->d_layers[0].vec_gridNum[1] << "\n"; 
for (int part_id =0; part_id < num_parts; part_id ++ )
{
file_stream_log << part_id << "\t";
file_stream_log <<this->d_layers[part_id].vec_size[0]<< "\t";
file_stream_log <<this->d_layers[part_id].vec_size[1]<< "\t";
file_stream_log <<this->d_layers[part_id].vec_gridNum[0]<< "\t";
file_stream_log <<this->d_layers[part_id].vec_gridNum[2]<< "\n";
}


file_stream_log.close();

 }
//generate the XDA file for the whole structure: 
// print the overall info: layer info
void GenerateMesh::print_info()
{
    std::cout << STRING_WRITE_BEGIN<< " print overall info" << endl;
    std::cout<< "+++ number of layers: " << d_no_layer << endl;
    std::cout<< "+++ for each layer" << endl;

    std::cout<< STRING_WRITE_END <<" print overall info" << endl;
} // print_info

