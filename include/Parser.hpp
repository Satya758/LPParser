
#ifndef LPPARSER_PARSER_HPP
#define LPPARSER_PARSER_HPP

#include "Problem.hpp"

#include "glpk.h"

namespace lp {

/**
 * Converts GLPK structures to Blaze matrices
 */
class Parser {
 private:
  glp_prob* _problem;

  glp_prob* createGlpProblem(const std::string filePath);

 public:
  Parser(const std::string filePath);

  ~Parser();

  Problem getProblem();
};
}  // lp

#endif  // LPPARSER_PARSER_HPP
