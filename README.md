# LP Parser
- LPParser is a C++ library to read IBM LP format files into [Blaze](https://bitbucket.org/blaze-lib/blaze) matrices or CCS (Compressed Column Storage) structures.
- LPParser relies on [GLPK](https://www.gnu.org/software/glpk/) to parse LP format files and converts GLPK structures to following form.
    ```
    min  c'*x
    s.t. A*x = b
         G*x <=_K h
    ```
 
# Dependencies
 - [Blaze](https://bitbucket.org/blaze-lib/blaze)
 - [GLPK](https://www.gnu.org/software/glpk/)
 - [SPDLOG](https://github.com/gabime/spdlog)
 - [CMake](cmake.org)
 - [Boost](boost.org), boost system and boost thread (blaze dependencies)
 - gcc 4.9.2
 
# Install
 - Generate make files passing 
   * -DGLPK_INC="" 
   * -DGLPK_LIB="" 
   * -DBLAZE_INC="" 
   * -DSPDLOG_INC="" 
   * -DBOOST_LIB=""
   * -DCMAKE_INSTALL_PREFIX=""
 - make install  
 
# Usage
 - Include header Parser.hpp and Problem.hpp
 - Link with parser static lib
 - Parser header has two methods which return problem data.
   ```
    getBlazeProblem()
    getCCSProblem()
   ```
 - Constructor of parser, first argument accepts lp format file and second argument if true converts all inequality constraints to equality constraints. (except implicit constraint x >= 0 and bounds).
 - If convertInequalityToEquality is true, number of rows in A will and increase in number of columns due to additional slack variables, but G matrix will be sparse (Only one column active for each row).
 - If convertInequalityToEquality is false, size of matrices are same as input as currently (Free and Fixed variables are not supported).
   
   ```
   Parser(const std::string filePath, const bool convertInequalityToEquality)
   ```
 - If you need to parse MPS files, adapt Parser.cpp constructor to read MPS files.   
  
# Limitations / TODO
 - Free and Fixed variables are not supported yet.
 - Cannot convert solution to original problem.
 - Only dealing with LP problems (No MILP or other cones).
 - Only tested in Ubuntu.
 
# Why / Motivation
 - There is no parser available to run new open source solvers like [ECOS](https://github.com/embotech/ecos), [SCS](https://github.com/cvxgrp/scs), [cvxopt](cvxopt.org) e.t.c.
 - Test your own implementation.

# Faq
 - Its not a presolver.