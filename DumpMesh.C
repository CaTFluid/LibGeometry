#include "InitialGeometry.h"
#include <string.h>
#include <iostream>
using namespace std;
using namespace LibGeometry;
using namespace MY_APP;

int main(int argc, char** argv)
{   
    string inputFile= "";
     if( argc < 2)   
    { 
     inputFile = "geometry.template";
     cout <<" default input:geometry.template " << endl;
     
    }
     else
    inputFile = argv[1];
    
    GenerateMesh my_mesh(inputFile);
    cout<< "end of individual part" <<endl;
    my_mesh.GenerateXDAFile();
    cout<< "end of whole structure" <<endl;
    return 0;
}
