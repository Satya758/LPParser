
#include "glpk.h"

#include "Parser.hpp"

namespace lp {

const std::shared_ptr<spdlog::logger> Parser::_logger =
    spdlog::stdout_logger_mt("console");

glp_prob* Parser::createGlpProblem(const std::string filePath) const {
  glp_prob* problem = glp_create_prob();

  glp_read_lp(problem, NULL, filePath.data());

  return problem;
}

Parser::Dimensions Parser::getFinalDimensions(int definedRows,
                                              int definedColumns) const {
  int equalityRows = 0;
  int inequalityRows = 0;
  int columns = definedColumns;

  size_t equalityNnz = 0;
  size_t inequalityNnz = 0;

  for (int r = 1; r <= definedRows; ++r) {
    int rowType = glp_get_row_type(_problem, r);
    int nnzRow = glp_get_mat_row(_problem, r, nullptr, nullptr);

    switch (rowType) {
      case GLP_FR:
        throw std::invalid_argument("There will be no free rows");
      case GLP_FX:
        ++equalityRows;
        equalityNnz += nnzRow;
        break;
      case GLP_UP:
      case GLP_LO:
        ++inequalityRows;
        inequalityNnz += nnzRow;
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
        // Positive orthant constraints
        ++inequalityRows;
        ++inequalityNnz;
        break;
      case GLP_DB:
        // What is format for double bound eg. low <= var <= high, or all
        // combinations!
        // eg. bin1 <= 300 is treated as double bound like 0 <= bin1 <= 300
        inequalityRows += 2;
        inequalityNnz += 2;
        break;
      default:
        throw std::invalid_argument("Wrong row type");
    }
  }

  _logger->info("Equality rows: {}, Inequality rows: {}, Columns: {}, "
                "Equality nnz: {}, Inequality nnz: {}",
                equalityRows, inequalityRows, columns, equalityNnz,
                inequalityNnz);
  // TODO Increase nnz by small factor to avoid any memory issues! But whats
  // the factor!?
  // TODO I think I have calculated nnz properly, but who knows
  return std::make_tuple(equalityRows, inequalityRows, columns, equalityNnz,
                         inequalityNnz);
}

Problem Parser::createProblem(Parser::Dimensions finalRowsAndColumns,
                              int definedRows, int definedColumns) const {
  Problem problem(std::get<0>(finalRowsAndColumns),
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

  A.reserve(std::get<3>(finalRowsAndColumns));
  G.reserve(std::get<4>(finalRowsAndColumns));

  _logger->info("Data reserved");

  for (int i = 1; i <= definedRows; ++i) {
    int nnzRow = glp_get_mat_row(_problem, i, index.get(), values.get());
    int rowType = glp_get_row_type(_problem, i);

    int sign = 1; // Default sign
    if (rowType == GLP_LO) {
      sign = -1;
    }

    for (int j = 1; j <= nnzRow; ++j) {

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

  problem.A = A;
  problem.G = G;

  return problem;
}

inline void Parser::createInequalityRowWithValueOne(
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

Problem Parser::getProblem() const {
  int definedRows = glp_get_num_rows(_problem);
  int definedColumns = glp_get_num_cols(_problem);

  Dimensions finalRowsAndColumns =
      getFinalDimensions(definedRows, definedColumns);

  return createProblem(finalRowsAndColumns, definedRows, definedColumns);
}

ClassicalProblem Parser::getClassicalProblem() const {
  _logger->info("Creation of Blaze objects started");
  const Problem problem = getProblem();
  _logger->info("Creation of Blaze objects Ended");

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

template <typename CP, typename RP, typename VAL>
void Parser::createSparseMatrix(CP& cp, RP& rp, VAL& value,
                                const SparseMatrix& coeffMatrix) const {
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

      ++columnPtr;

      value[rowPtr] = iter.value();
      rp[rowPtr] = iter.index();

      ++rowPtr;
    });

    cp[i + 1] = columnPtr;
  }
}

} // lp