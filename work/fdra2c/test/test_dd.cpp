#include <iostream>
#include "src/DDouble.hpp"
using namespace fdra;
using namespace std;

int main(int argc, char** argv)
{
  DDouble t = 1.0;
  double h = 1.0e-19;
  DDouble tnew = t + h;
  h = tnew - t;
  cout << "h = " << h << endl;
  h = t - tnew;
  cout << "h = " << h << endl;
  cout << "1.1 - t = " << (1.1 - t).ToDouble() << endl;
}
