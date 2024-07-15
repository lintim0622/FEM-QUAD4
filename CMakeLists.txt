cmake_minimum_required (VERSION 3.8)

project ("FEM project")

include_directories("C://Program Files")
# ${CMAKE_CURRENT_SOURCE_DIR}/googletest/include

# list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/googletest)

# 將來源新增到此專案的可執行檔。
add_executable (CMakeTarget "mesh.cpp" 
                            "mesh.h" 
                            "main.cpp"
                            "unit_test.cpp" 
                            "shape_function.h"
                            "shape_function.cpp" 
                            "calculate.h"
                            "calculate.cpp"
                             "static_analysis.h" "static_analysis.cpp" "read_initial.h" "reda_initial.cpp")
                            #  
                            
                            

# target_link_libraries(CMakeTarget gtest)

# 使用 FetchContent 模組
include(FetchContent)
# 下載和引入 Google Test
# FetchContent_Declare(
#     googletest
#     GIT_REPOSITORY https://github.com/google/googletest.git
#     GIT_TAG release-1.11.0
# )
# FetchContent_MakeAvailable(googletest)