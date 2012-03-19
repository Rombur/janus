#include <cassert>
#include <string>
#include <TRIANGULATION.hh>

using namespace std;

int main(int argc,char** argv)
{
  assert(argc==4);

  string geometry_filename(argv[1]);
  string parameters_filename(argv[2]);
  string output_filename(argv[3]);

  TRIANGULATION triangulation(&geometry_filename);

  return 0;
}
