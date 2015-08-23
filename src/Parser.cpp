

#include "Parser.hpp"

glp_prob* lp::Parser::createGlpProblem(const std::string filePath) {
  glp_prob* problem = glp_create_prob();

  glp_read_lp(problem, NULL, filePath.data());

  return problem;
}

lp::Parser::Parser(const std::string filePath)
    : _problem(createGlpProblem(filePath)) {}

lp::Problem lp::Parser::getProblem() {}