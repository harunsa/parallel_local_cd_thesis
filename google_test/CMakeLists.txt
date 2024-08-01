cmake_minimum_required(VERSION 3.26)

# 'Google_test' is the subproject name
project(google_test)

# 'lib' is the folder with Google Test sources
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# 'Google_Tests_run' is the target name
# 'test1.cpp test2.cpp' are source files with tests
add_executable(Google_Test_run test1.cpp)
target_link_libraries(Google_Test_run gtest gtest_main)

target_link_libraries(Google_Test_run
        adj_list
)

target_include_directories(Google_Test_run PUBLIC
        ../include/adj_list
)