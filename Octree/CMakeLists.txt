cmake_minimum_required(VERSION 3.15)
project(cmake_vtk)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(VTK COMPONENTS
  CommonColor
  CommonCore
  FiltersSources
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

include_directories($(CMAKE_CURRENT_SOURCE_DIR)/source)

add_executable(cmake_vtk_main source/main.cpp)

target_link_libraries(cmake_vtk_main PRIVATE ${VTK_LIBRARIES} )


# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS cmake_vtk_main cmake_vtk_readerfile cmake_vtk_drawcube
  MODULES ${VTK_LIBRARIES}
)

