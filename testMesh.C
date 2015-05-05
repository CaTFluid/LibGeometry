/***********************************************************************
	EXEMPLO 1.2
	This example introduces Elem and Node classes
************************************************************************/

#include "equation_systems.h"
#include "numeric_vector.h"
#include <vector>
#include <string>
#include <sstream> 
// C++ headers
#include <iostream>

// libMesh headers
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"

#include "elem.h"
#include "node.h"
#include <stdio.h>
#include "exodusII_io.h"
#include "gmsh_io.h"
using namespace std;
int main (int argc, char** argv)
{
	// Initialize the library (obligatory instruction for every programa).
	LibMeshInit init(argc, argv);
	int nDim = 3;
	string meshFile = "exA.xda";
	if (argc == 1)
	{
	 nDim =3;
	 string meshFile = "exA.xda";
	 cout<<"default: 3D mesh" << endl;	
	 cout<< "default: input xdafile: " << meshFile << endl;
	}
	else if (argc ==2)
	{
	 nDim =3;
	 cout <<"default: 3D mesh" << endl;
	 meshFile = argv[1];
	 cout << "input xdafile:" << meshFile << endl;
	}
	else if (argc==3)
	{
	 nDim = atoi (argv[2]);
	 cout <<"default: 3D mesh" << endl;
	 meshFile = argv[1];
	 cout << "input xdafile:" << meshFile << endl;
	}
	
	//int nDim =3;
	// Create a Mesh object for a 2d domain.
	Mesh mesh(nDim);
	// READ the mesh file
	//string meshFile = "exA.xda";
	// Print out info of the mesh file
	mesh.read(meshFile);
	// OutPut the mesh file for visualization
	mesh.print_info();
        //mesh.write("a.vtu");
	string StrMSH = meshFile.substr(0, meshFile.size() - 4) + ".msh";
	mesh.write(StrMSH);


	// Iterator used to percorrer xxx the elements of the mesh, similar to the iterators of the C++ containers (set, list, etc)
	Mesh::element_iterator       it_el      = mesh.active_local_elements_begin();//mesh.elements_begin();
	const Mesh::element_iterator it_last_el = mesh.active_local_elements_end();//mesh.elements_end();
	for ( ; it_el != it_last_el ; ++it_el) {

		// The iterator gives us a pointer to the element (not the element object in itself); We get the pointer from the iterator only to make the syntex more clear;
		const Elem* elem = *it_el;

		// Let's show some information about this element (note: the numbering in libmesh data structures always start in 0, so for the output we add 1)
		cout << "Element number: " << elem->id()+1 << endl;
		cout << "Element subdomiain id: " << elem->subdomain_id () << endl;

		// percorrer the element nodes
		for(unsigned int i=0; i < elem->n_nodes(); i++)  {
			// Show some informations about the node whose local number is "i"; 
			//NOTE: the "get_node" member of the classe "Elem" gives us a pointer to the Node object, so we can do everything we want with that node (for example, change it's coordinates, )
			Node* no_temp = elem->get_node(i);
			cout << "	node local number: " << i+1 << "	node global number: " << no_temp->id() + 1 << "	node coordinates: " << *no_temp;
		}

		// find who the neighbors of this element are
		cout << endl;
		for(unsigned int i=0; i < elem->n_sides(); i++)  {
			const Elem* elem_neighbor = elem->neighbor(i);

			//if this side is on the boundary, we get a NULL pointer
			if(elem_neighbor == NULL)
				cout << "	side " << i+1 << " does not have a neighbor (side on border)" << endl;
			else
				cout << "	side " << i+1 << " has element number " << elem_neighbor->id()+1 << " as neighbor" << endl;
		}
		cout << endl;
		// added 02-06 --> write subdomain id to see how it outputs:
		Elem* non_constant_elem = *it_el;
		libMesh::subdomain_id_type  & id1= non_constant_elem->subdomain_id();
		cout << "subdomain_id before modify: " << id1 <<endl;
		if (non_constant_elem->id() + 1 > 2)
			{ 
			id1 = id1+1;
			cout << "subdomain_id after modify: " << non_constant_elem->subdomain_id() << endl;
			}
			 //id1=1;}
		 else	{ id1 = id1 + 2;
			 cout << "subdomain_id after modify: " << non_constant_elem->subdomain_id() << endl;
			}
		
		
	}
}
