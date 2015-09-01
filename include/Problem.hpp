#ifndef PROBLEM_HPP
#define PROBLEM_HPP

/**
 * Subset of file copied from LinearProgramming repository Problem definiton
 */
#include <memory>
#include <string>

#include <blaze/Blaze.h>

// TODO namespace lp is used as working name until I find name I like :-)
namespace lp {

// If double is not good enough for all platforms take it as template parameter
using SparseMatrix = blaze::CompressedMatrix<double, blaze::columnMajor>;

// FIXME If not used remove this definition, as for all RHS vectors dense
// vectors are considered
using SparseVector = blaze::CompressedVector<double, blaze::columnVector>;

// All vectors used are dense vectors
using DenseVector = blaze::DynamicVector<double, blaze::columnVector>;

/**
 *Defines the problem and options to solver
 *
 *Solves primal problem
 * 	minimize   c'x
 *    subject to   Gx + s = h or Gx <= h (Inequality constraints)
 * 		   Given Inequality constraints should be of the form Gx <= h
 * 		   (s is added to simplify algorithm)
 * 		   Ax = b (Equality constraints)
 *		   s >= 0, s is subset of Cone
 * 	Dimensions c = n x 1
 * 		   x = n x 1
 * 		   G = m x n
 *		   h = m x 1 (Slack variable s has same dimension as h)
 *		   s = m x 1 (Slack variable is not passed from user)
 *		   A = me x n (me is m rows for equality constraints)
 * 		   b = me x 1
 *
 *And its dual problem
 * 	maximize    -h'z - b'y
 * 	subject to  G'z + A'y + c = 0
 *		    z >= 0, z is subset of Cone
 *
 *Cone is cartesian product of multiple cones ( Positive orthant, Second order
 *cone, SDP)
 *
 */
class Problem {
public:
  // Options are defaulted
  Problem(size_t equalityRows, size_t inequalityRows, size_t columns)
      : equalityRows(equalityRows),
        inequalityRows(inequalityRows),
        columns(columns),
        c(this->columns),
        G(this->inequalityRows, this->columns),
        h(this->inequalityRows),
        A(this->equalityRows, this->columns),
        b(this->equalityRows) {}

  const size_t equalityRows;
  const size_t inequalityRows;
  const size_t columns;

  DenseVector c; // Objective coefficients

  SparseMatrix G; // Inequality constraints
  DenseVector h;  // RHS of Inequality constraints

  SparseMatrix A; // Equality constraints
  DenseVector b;  // RHS of equality constraints

  friend std::ostream& operator<<(std::ostream& os, const Problem& problem);
};

inline std::ostream& operator<<(std::ostream& os, const Problem& problem) {
  std::cout << "Equality matrix: \n" << problem.A << std::endl;
  std::cout << "Inequality matrix: \n" << problem.G << std::endl;

  std::cout << "Equality vector: \n" << problem.b << std::endl;
  std::cout << "Inequality vector: \n" << problem.h << std::endl;

  std::cout << "Objective : \n" << problem.c << std::endl;

  return os;
}

template <class T> using Array = std::unique_ptr<T[]>;

/**
 * Plain pointers instead of blaze objects.
 *
 * not first class citizen, created from blaze.
 */
class ClassicalProblem {
public:
  ClassicalProblem(const size_t equalityRows, const size_t inequalityRows,
                   const size_t columns, const size_t Gnnz, const size_t Annz)
      : equalityRows(equalityRows),
        inequalityRows(inequalityRows),
        columns(columns),
        c(std::make_unique<double[]>(columns)),
        GColumnptr(std::make_unique<size_t[]>(columns + 1)),
        GRowPtr(std::make_unique<size_t[]>(Gnnz)),
        GValue(std::make_unique<double[]>(Gnnz)),
        AColumnptr(std::make_unique<size_t[]>(columns + 1)),
        ARowPtr(std::make_unique<size_t[]>(Annz)),
        AValue(std::make_unique<double[]>(Annz)),
        h(std::make_unique<double[]>(inequalityRows)),
        b(std::make_unique<double[]>(equalityRows)) {}

  const size_t equalityRows;
  const size_t inequalityRows;
  const size_t columns;

  Array<double> c;

  Array<size_t> GColumnptr;
  Array<size_t> GRowPtr;
  Array<double> GValue;

  Array<size_t> AColumnptr;
  Array<size_t> ARowPtr;
  Array<double> AValue;

  Array<double> h;
  Array<double> b;

  friend std::ostream& operator<<(std::ostream& os,
                                  const ClassicalProblem& problem);
};

template <typename T>
void printArray(const std::string& context, size_t size,
                const Array<T>& array) {
  std::cout << "\n\n" << context << std::endl;

  for (size_t j = 0; j < size; ++j) {
    std::cout << array[j] << ", ";
  }
}

inline std::ostream& operator<<(std::ostream& os,
                                const ClassicalProblem& problem) {
  printArray("Objective matrix", problem.columns, problem.c);

  printArray("Inequality matrix Column ptr", problem.columns + 1,
             problem.GColumnptr);
  // Get nnz from column ptr
  printArray("Inequality matrix Row ptr", problem.GColumnptr[problem.columns],
             problem.GRowPtr);
  printArray("Inequality matrix values", problem.GColumnptr[problem.columns],
             problem.GValue);

  printArray("Equality matrix Column ptr", problem.columns + 1,
             problem.AColumnptr);
  // Get nnz from column ptr
  printArray("Equality matrix Row ptr", problem.AColumnptr[problem.columns],
             problem.ARowPtr);
  printArray("Equality matrix values", problem.AColumnptr[problem.columns],
             problem.AValue);

  printArray("Equality RHS", problem.equalityRows, problem.b);

  printArray("Inequality RHS", problem.inequalityRows, problem.h);

  return os;
}
} // lp
#endif // PROBLEM_HPP