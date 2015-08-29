
#ifndef LPPARSER_PARSER_HPP
#define LPPARSER_PARSER_HPP

#include <memory>
#include <utility>
#include <algorithm>
#include <iostream>

#include "Problem.hpp"

#include "glpk.h"

#include "spdlog/spdlog.h"

namespace lp {

/**
 * Converts GLPK structures to Blaze matrices
 *
 * FIXME use size_t
 */
class Parser {
private:
  const static std::shared_ptr<spdlog::logger> _logger;

  using Dimensions = std::tuple<int, int, int>;

  glp_prob* _problem;

  bool _convertInequalityToEquality = false;

  glp_prob* createGlpProblem(const std::string filePath) {
    glp_prob* problem = glp_create_prob();

    glp_read_lp(problem, NULL, filePath.data());

    return problem;
  }

  /**
   * Returns final Equality rows and Inequality rows, used to initialize
   * SparseMatrix and DenseVector sizes
   */
  Dimensions getFinalRowsAndColumns(int definedRows, int definedColumns) {
    int equalityRows = 0;
    int inequalityRows = 0;
    int columns = definedColumns;

    for (int r = 1; r <= definedRows; ++r) {
      int rowType = glp_get_row_type(_problem, r);

      // _logger->info("Row type {}", rowType);

      switch (rowType) {
        case GLP_FR:
          throw std::invalid_argument("There will be no free rows");
        case GLP_FX:
          ++equalityRows;
          break;
        case GLP_UP:
        case GLP_LO:
          ++inequalityRows;
          break;
        case GLP_DB:
          throw std::invalid_argument(
              "Double bounds in rows, lazy to implement as they are rare");
        default:
          throw std::invalid_argument("Wrong row type");
      }
    }

    for (int c = 1; c <= definedColumns; ++c) {
      int columnType = glp_get_col_type(_problem, c);

      //      _logger->info("Column type {}", columnType);

      switch (columnType) {
        case GLP_FR:
          // Free variable results in two additional variables eg. x = u - v (x
          // is free variable and u, v are any positive variables ), this
          // increases columns and inequality rows by two
          // inequalityRows += 2; // TODO Deal later
          // Column is replaced by two columns so increase column count by one
          // ++columns; // TODO Deal later
          throw std::invalid_argument("Free variables are dealt later");
          break;
        case GLP_FX:
          throw std::invalid_argument(
              "Fixed variables are constants are not considered now, until I "
              "face a problem :-)");
        case GLP_UP:
        case GLP_LO:
          ++inequalityRows;
          break;
        case GLP_DB:
          // What is format for double bound eg. low <= var <= high, or all
          // combinations!
          // eg. bin1 <= 300 is treated as double bound like 0 <= bin1 <= 300
          inequalityRows += 2;
          break;
        default:
          throw std::invalid_argument("Wrong row type");
      }
    }

    return std::make_tuple(equalityRows, inequalityRows, columns);
  }

  lp::Problem createProblem(Dimensions finalRowsAndColumns, int definedRows,
                            int definedColumns) {
    lp::Problem problem(std::get<0>(finalRowsAndColumns),
                        std::get<1>(finalRowsAndColumns),
                        std::get<2>(finalRowsAndColumns));

    // Used to increment equality matrix rows, keeps track of equalityRows
    // inserted
    int equalityRow = 0;
    int inequalityRow = 0;

    // As glp index starts with 1 we have added +1 memory
    std::unique_ptr<int[]> index;
    index.reset(new int[definedColumns + 1]);

    std::unique_ptr<double[]> values;
    values.reset(new double[definedColumns + 1]);

    // Creating new row matrix not using Problems object as they are column
    // major
    blaze::CompressedMatrix<double> A(std::get<0>(finalRowsAndColumns),
                                      std::get<2>(finalRowsAndColumns));
    blaze::CompressedMatrix<double> G(std::get<1>(finalRowsAndColumns),
                                      std::get<2>(finalRowsAndColumns));

    problem.b.reserve(problem.b.size());
    problem.h.reserve(problem.h.size());

    for (int i = 1; i <= definedRows; ++i) {
      int nnzRow = glp_get_mat_row(_problem, i, index.get(), values.get());
      int rowType = glp_get_row_type(_problem, i);

      if (rowType == GLP_FX) {
        A.reserve(equalityRow, nnzRow);
      } else if (rowType == GLP_UP || rowType == GLP_LO) {
        G.reserve(inequalityRow, nnzRow);
      }

      int sign = 1; // Default sign
      if (rowType == GLP_LO) {
        sign = -1;
      }

      for (int j = 1; j <= nnzRow; ++j) {
        /*_logger->info("Row index: {}, Column index: {}, Column value: {}", i,
                      index[j], values[j]);*/

        if (rowType == GLP_FX) {
          A.append(equalityRow, index[j] - 1, values[j]);
          problem.b[equalityRow] = glp_get_row_ub(_problem, i);
        } else if (rowType == GLP_UP || rowType == GLP_LO) {
          G.append(inequalityRow, index[j] - 1, sign * values[j]);

          // Either rowType or sign computed above can be used
          if (rowType == GLP_UP) {
            problem.h[inequalityRow] = glp_get_row_ub(_problem, i);
          } else if (rowType == GLP_LO) {
            problem.h[inequalityRow] = -glp_get_row_lb(_problem, i);
          }
        }
      }

      if (rowType == GLP_FX) {
        A.finalize(equalityRow++);
      } else if (rowType == GLP_UP || rowType == GLP_LO) {
        G.finalize(inequalityRow++);
      }
    }

    for (int columnIndex = 1; columnIndex <= definedColumns; ++columnIndex) {
      int columnType = glp_get_col_type(_problem, columnIndex);

      if (columnType == GLP_FR) {
        throw std::invalid_argument("Free variables are dealt later");
      } else if (columnType == GLP_UP) {
        // This wont happen, if there is anything it would turn up in GLP_DB
      } else if (columnType == GLP_LO) { // These are x >= 0, default bound
        createInequalityRowWithValueOne(G, inequalityRow, columnIndex - 1, -1,
                                        problem.h); // Lower bound
      } else if (columnType == GLP_DB) {
        createInequalityRowWithValueOne(G, inequalityRow, columnIndex - 1, 1,
                                        problem.h); // Upper bound
        createInequalityRowWithValueOne(G, inequalityRow, columnIndex - 1, -1,
                                        problem.h); // Lower bound
      }

      // Objective function
      problem.c[columnIndex - 1] = glp_get_obj_coef(_problem, columnIndex);
    }

    // FIXME _logger Not wokring why?
    //    _logger->info("Equality rows: {}", A);
    //    _logger->info("Ineqaulity rows: {}", G);

    //    _logger->info("Inequality rows: ") << G;

    problem.A = std::move(A);
    problem.G = std::move(G);

    return problem;
  }

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
      blaze::DynamicVector<double, blaze::columnVector>& h) const {
    if (value == -1) { // Lower bound
      h[row] = -glp_get_col_lb(_problem, columnIndex + 1);
    } else if (value == 1) { // Upper bound
      h[row] = glp_get_col_ub(_problem, columnIndex + 1);
    }

    G.reserve(row, 1);
    G.append(row, columnIndex, value);
    G.finalize(row++);
  }

public:
  Parser(const std::string filePath, const bool convertInequalityToEquality)
      : _problem(createGlpProblem(filePath)),
        _convertInequalityToEquality(convertInequalityToEquality) {}

  ~Parser() { glp_free(_problem); };

  Problem getProblem() {
    _logger->info("Parsing started");

    int definedRows = glp_get_num_rows(_problem);
    int definedColumns = glp_get_num_cols(_problem);

    Dimensions finalRowsAndColumns =
        getFinalRowsAndColumns(definedRows, definedColumns);

    _logger->info("Equality rows {}, Inequality rows {}, Columns {}",
                  std::get<0>(finalRowsAndColumns),
                  std::get<1>(finalRowsAndColumns),
                  std::get<2>(finalRowsAndColumns));

    return createProblem(finalRowsAndColumns, definedRows, definedColumns);
  }

  ClassicalProblem getClassicalProblem() {
    const Problem problem = getProblem();

    ClassicalProblem cProblem(problem.equalityRows, problem.inequalityRows,
                              problem.columns, problem.G.nonZeros(),
                              problem.A.nonZeros());
    // Objective
    for (size_t j = 0; j < problem.c.size(); ++j) {
      cProblem.c[j] = problem.c[j];
    }

    // RHS Equality
    for (size_t k = 0; k < problem.b.size(); ++k) {
      cProblem.bValue[k] = problem.b[k];
    }
    // RHS Inequality
    for (size_t k = 0; k < problem.h.size(); ++k) {
      cProblem.hValue[k] = problem.h[k];
    }

    // G
    createSparseMatrix(cProblem.GColumnptr, cProblem.GRowPtr, cProblem.GValue,
                       problem.G);
    // A
    createSparseMatrix(cProblem.AColumnptr, cProblem.ARowPtr, cProblem.AValue,
                       problem.A);

    return cProblem;
  }

  /**
   * Creates regular SparseMatrix from blaze SparseMatrix
   */
  template <typename CP, typename RP, typename VAL>
  void createSparseMatrix(CP& cp, RP& rp, VAL& value,
                          const SparseMatrix& coeffMatrix) {
    // FIXME lambda function does not allow pointer as input, may be problem
    // with blaze iterator
    using ColumnIterator =
        typename std::remove_pointer<SparseMatrix::ConstIterator>::type;

    // Nature of column pointer, always starts with 0 and ends with nnz
    cp[0] = 0;
    int columnPtr = 0;
    int rowPtr = 0;

    for (size_t i = 0; i < coeffMatrix.columns(); ++i) {
      std::for_each(blaze::cbegin(coeffMatrix, i), blaze::cend(coeffMatrix, i),
                    [&](const ColumnIterator& iter) {
        /*_logger->info("Column: {}, Row: {}, Value: {}", i, iter.index(),
                      iter.value());*/

        ++columnPtr;

        value[rowPtr] = iter.value();
        rp[rowPtr] = iter.index();

        /*_logger->info("Inserted value: {} at rowPtr {}", value[rowPtr],
         * rowPtr);*/

        ++rowPtr;
      });

      cp[i + 1] = columnPtr;
    }
  }
};

const std::shared_ptr<spdlog::logger> Parser::_logger =
    spdlog::stdout_logger_mt("console");

} // lp

#endif // LPPARSER_PARSER_HPP
