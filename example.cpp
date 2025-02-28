// define nil (how to circumvent that define?)
#define nil NULL
// necessary NEURON includes

#include "./nrn/src/ivoc/ivocmain.cpp"
#include "./nrn/src/ivoc/ocjump.cpp"
// necessary stdlib includes
#include <string>
#include <vector>
#include <sstream>
#include <map>

// extern C functions
extern bool hoc_valid_stmt(const char* stmt, Object* ob);
extern int ivocmain(int, char**, char**);
extern double hoc_ac_;

// oc environment
static char* env[] = {0};



// main
int main(int argc, char** argv) {
   const int init = 0;
   ivocmain(init, argv, env);
   const char* stmt = "x = 5";
   const char* stmt2 = "hoc_ac_ = x";
   const bool ret = hoc_valid_stmt(stmt, 0);
   return 1;
}
