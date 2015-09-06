
#ifndef LPPARSER_PARSER_HPP
#define LPPARSER_PARSER_HPP

#include <memory>
#include <utility>
#include <algorithm>
#include <iostream>

#include "Problem.hpp"

#include "glpk.h"

#include "spdlog/spdlog.h"

namespace lpp {

/**
 * Converts GLPK structures to Blaze matrices
 *
 * FIXME use size_t
 */
class Parser {
private:
  const static std::shared_ptr<spdlog::logger> _logger;

  using Dimensions = std::tuple<int, int, int, int, int>;

  glp_prob* _problem;

  // flag to indicate to convert Inequality constraints to equality constraints
  // by adding slack variable, so there will additional column for each
  // inequality constraint, sign of slack depends on inequality operator
  bool _convertICToEC = false;

  glp_prob* createGlpProblem(const std::string filePath) const;

  /**
   * Returns final Equality rows and Inequality rows, used to initialize
   * SparseMatrix and DenseVector sizes
   */
  Dimensions getFinalDimensions(int definedRows, int definedColumns) const;

  Problem createProblem(Dimensions finalRowsAndColumns, int definedRows,
                        int definedColumns) const;

  /**
   * value though double, all we need in this case is 1, so int. value is used
   * as sign
   *
   * value 1 is upper bound
   * value -1 is lower bound
   */
  void createInequalityRowWithValueOne(
      blaze::CompressedMatrix<double>& G, int& row, const int columnIndex,
      const int value,
      blaze::DynamicVector<double, blaze::columnVector>& h) const;

  /**
   * Creates regular SparseMatrix from blaze SparseMatrix
   */
  template <typename CP, typename RP, typename VAL>
  void createSparseMatrix(CP& cp, RP& rp, VAL& value,
                          const SparseMatrix& coeffMatrix) const;

public:
  Parser(const std::string filePath, const bool convertInequalityToEquality)
      : _problem(createGlpProblem(filePath)),
        _convertICToEC(convertInequalityToEquality) {}

  ~Parser() { glp_free(_problem); };

  /**
   *
   */
  Problem getBlazeProblem() const;

  /**
   * get Column compressed storage
   */
  CCSProblem getCCSProblem() const;
};

} // lp

#endif // LPPARSER_PARSER_HPP
