# CmakeLists.txt in Tools

find_package(Qt6 REQUIRED COMPONENTS Core Widgets)
set(QT_USE_QTSVG TRUE)

file(GLOB Tools_CPP "*.cpp")
file(GLOB Tools_H "*.h")
file(GLOB Tools_UI "*.ui")

qt6_wrap_ui(Tools_H ${Tools_UI})

list(APPEND Tools_Sources ${Tools_CPP} ${Tools_H})

# This is for moc files
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(CornucopiaTools ${Tools_Sources})

set_target_properties(CornucopiaTools PROPERTIES AUTOMOC TRUE)

# And this is for config.h
include_directories(${CMAKE_BINARY_DIR})

include_directories(${Cornucopia_SOURCE_DIR}/Cornucopia)

target_link_libraries(CornucopiaTools PRIVATE Cornucopia Qt6::Core Qt6::Widgets)

configure_file(config.h.cmake config.h)
