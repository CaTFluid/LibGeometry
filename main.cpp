#include "InitialGeometry.h"
#include <string.h>
#include <iostream>
using namespace std;
using namespace LibGeometry;
using namespace MY_APP;

int main()
{   string inputFile= "geometry.template";
    GenerateMesh my_mesh(inputFile);
    cout<< "end" <<endl;

    return 1;
}
