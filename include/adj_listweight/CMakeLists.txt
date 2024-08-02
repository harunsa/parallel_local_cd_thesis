cmake_minimum_required(VERSION 3.26)
project(adj_listweight)

add_library(adj_listweight STATIC
        adj_listweight.cpp
)

set_target_properties(adj_listweight PROPERTIES PUBLIC_HEADER adj_listweight.h)

target_link_libraries(adj_listweight PUBLIC
        libcuckoo
        parlay
        bulk
)

target_include_directories(adj_listweight PUBLIC
        ../libcuckoo/libcuckoo
        ../parlaylib/
        ../ext/bulk
)