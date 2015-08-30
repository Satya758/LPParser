cmake_minimum_required(VERSION 3.3)
project(LPParser)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -std=c++14 -Wall")

add_subdirectory(include)

include_directories(include)

# GLPK, Only parser related functions are used
set(GLPK_INC "" CACHE PATH "GLPK library, only parser functions are used")
set(GLPK_LIB "" CACHE PATH "GLPK library, only parser functions are used")

# Blaze, core linear algebra
set(BLAZE_INC "" CACHE PATH "Blaze header library")
# Logging
set(SPDLOG_INC "" CACHE PATH "Logging library")
# Boost, dependent for blaze
set(BOOST_LIB "" CACHE PATH "Boost library")
find_library(BOOST_SYSTEM boost_system PATHS ${BOOST_LIB} NO_DEFAULT_PATH)
find_library(BOOST_THREAD boost_thread PATHS ${BOOST_LIB} NO_DEFAULT_PATH)

# Name of lib and Variable being same does not work so added GLPK_P instead of GLPK
find_library(GLPK_P glpk PATHS ${GLPK_LIB})
include_directories(SYSTEM ${GLPK_INC} SYSTEM ${BLAZE_INC} SYSTEM ${SPDLOG_INC})

file(GLOB parser_src src/*.cpp)
add_library(parser ${parser_src})

# TODO Used for testing
set(SOURCE_FILES main.cpp)
add_executable(LPParser ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME}
        ${BOOST_SYSTEM} ${BOOST_THREAD}
        parser
        ${GLPK_P}
        pthread)

install(TARGETS parser DESTINATION lib)
install(FILES include/Parser.hpp include/Problem.hpp DESTINATION include)
