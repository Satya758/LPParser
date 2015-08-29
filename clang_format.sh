#!/bin/bash

shopt -s globstar
shopt -s nullglob
shopt -s nocaseglob


for c_file in $(find ./ -path ./build -prune -o -name '*.cpp' -or -name '*.hpp')
do
  clang-format-3.5 -style=Mozilla -i "$c_file"
  echo "$c_file"
done
