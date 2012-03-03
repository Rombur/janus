#include <cassert>
#include <string>

using namespace std;

int main(int argc,char** argv)
{
  assert(argc==4);

  string output_filename(argv[3]);

  return 0;
}
