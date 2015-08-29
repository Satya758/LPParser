
#include "glpk.h"

#include "Parser.hpp"

/**
 * TODO Currently file is ignored for sake ease of develoing
 */

glp_prob* lp::Parser::createGlpProblem(const std::string filePath) {
  glp_prob* problem = glp_create_prob();

  glp_read_lp(problem, NULL, filePath.data());

  _logger->info("Dead in water: {}", glp_get_num_rows(problem));
  return problem;
}

lp::Parser::Parser(const std::string filePath)
    : _problem(createGlpProblem(filePath)),
      _logger(spdlog::stdout_logger_mt("console")) {}

lp::Problem lp::Parser::getProblem() {
  //   int rows = glp_get_num_rows(_problem);

  //  _logger->info("Number of rows: {}", rows);

  //  return nullptr;

  lp::Problem problem(5, 5, 5);

  return problem;
}

// lp::Parser::~Parser() { /*glp_free(_problem);*/ }

lp::DenseVector getObjective(const glp_prob* problem) {}