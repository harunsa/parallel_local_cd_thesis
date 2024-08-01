cmake_minimum_required(VERSION 3.26)
project(adj_list)

add_library(adj_list STATIC
        adj_list.cpp
)

set_target_properties(adj_list PROPERTIES PUBLIC_HEADER adj_list.h)

target_link_libraries(adj_list PUBLIC
        libcuckoo
        parlay
)

target_include_directories(adj_list PUBLIC
        ../libcuckoo/libcuckoo
        ../parlaylib/
)