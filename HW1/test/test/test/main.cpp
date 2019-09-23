

#include <iostream>
#include <fstream>
using namespace std;

int main( int argc, char* argv[] )
{
    ofstream myfile;
    myfile.open ("example.csv");
    myfile << 100 << " test "  << endl;
}
