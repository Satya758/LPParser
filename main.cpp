#include <iostream>

#include "Parser.hpp"

using namespace std;
int main() {
  lp::Parser parser("/home/satya/LPProblems/QiTest.lp", false);
  //  lp::Parser parser("/home/satya/LPProblems/test/afiro.lp", false);

  lp::Problem problem = parser.getProblem();
  //
  cout << problem << endl;

  lp::ClassicalProblem cProblem = parser.getClassicalProblem();

  cout << cProblem << endl;

  return 0;
}