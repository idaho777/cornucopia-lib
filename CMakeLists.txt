# Root Cornucopia CMakeLists.txt

cmake_minimum_required(VERSION 3.10)

project(Cornucopia)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build, options are: \
        None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if(MSVC)
  add_definitions(/arch:SSE2 /fp:fast)
  # else() add_definitions(-ffast-math)
endif()

# Find Eigen 3
set(CMAKE_PREFIX_PATH ${Cornucopia_SOURCE_DIR}/../ ${CMAKE_PREFIX_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${Cornucopia_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(Cornucopia)
add_subdirectory(DemoUI)
add_subdirectory(Tools)
add_subdirectory(Test)

include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_NAME, "Cornucopia")
set(CPACK_PACKAGE_VENDOR "Ilya Baran")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/gpl.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_EXECUTABLES "DemoUI" "Cornucopia DemoUI")
include(CPack)
