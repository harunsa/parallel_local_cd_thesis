cmake_minimum_required(VERSION 3.26)
project(tempus)

set(CMAKE_CXX_STANDARD 17)



add_executable(tempus main.cpp)

add_subdirectory(include/libcuckoo)
add_subdirectory(include/parlaylib/)
add_subdirectory(google_test)
add_subdirectory(include/dyncomm)
add_subdirectory(include/adj_listweight)
add_subdirectory(include/eigen)
add_subdirectory(include/spectra)
add_subdirectory(ext/bulk)


target_link_libraries(tempus PUBLIC
        adj_listweight
        eigen
        dyncomm
        bulk_thread
)

target_include_directories(tempus PUBLIC
        include/adj_listweight
        include/eigen
        include/dyncomm
        include/bulk
)

