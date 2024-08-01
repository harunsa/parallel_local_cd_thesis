cmake_minimum_required(VERSION 3.26)
project(dyncomm)

add_library(dyncomm STATIC
        DynCommDriver.cpp
        RealPruneContainer.cpp
        PruneRunner.cpp
        Prune.cpp
        TemporalWeightedMinhashGadget.cpp
        TemporalNode.cpp
)

set_target_properties(dyncomm PROPERTIES PUBLIC_HEADER dyncomm.h)

target_link_libraries(dyncomm PUBLIC
        libcuckoo
        parlay
        eigen
        adj_listweight
        bulk
#       spectra
)

target_include_directories(dyncomm PUBLIC
        ../libcuckoo/libcuckoo
        ../parlaylib/
        ../eigen/Eigen
        ../adj_listweight
        ../spectra
        ../bulk
)