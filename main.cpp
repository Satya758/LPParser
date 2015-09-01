#include <iostream>

#include "Parser.hpp"
#include "spdlog/spdlog.h"

using namespace std;
int main() {

  const std::shared_ptr<spdlog::logger> _logger =
      spdlog::stdout_logger_mt("mainConsole");

  _logger->info("Parsing started");

//  lp::Parser parser("/home/satya/LPProblems/QiTest.lp", true);
//    lp::Parser parser("/home/satya/LPProblems/test/afiro.lp", true);
//    lp::Parser parser("/home/satya/LPProblems/test/maros-r7.lp", false);
  // TODO Free and Fixed variables are not handled yet
  //  lp::Parser parser("/home/satya/LPProblems/test/test.lp", false);

  lp::Parser parser("/home/satya/LPProblems/test/blend.lp", true);

  lp::Problem problem = parser.getProblem();
  //
  cout << problem << endl;

  //  lp::ClassicalProblem cProblem = parser.getClassicalProblem();

  _logger->info("Parsing ended");

  //    cout << cProblem << endl;

  return 0;
}