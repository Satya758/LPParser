#ifndef PROBLEM_HPP
#define PROBLEM_HPP

/**
 * Subset of file copied from LinearProgramming repository Problem definiton
 */

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
  Problem(int equalityRows, int inequalityRows, int columns)
      : equalityRows(equalityRows),
        inequalityRows(inequalityRows),
        columns(columns),
        c(this->columns),
        G(this->inequalityRows, this->columns),
        h(this->inequalityRows),
        A(this->equalityRows, this->columns),
        b(this->equalityRows) {}

  const int equalityRows;
  const int inequalityRows;
  const int columns;

  DenseVector c;  // Objective coefficients

  SparseMatrix G;  // Inequality constraints
  DenseVector h;   // RHS of Inequality constraints

  SparseMatrix A;  // Equality constraints
  DenseVector b;   // RHS of equality constraints
};

}  // lp
#endif  // PROBLEM_HPP