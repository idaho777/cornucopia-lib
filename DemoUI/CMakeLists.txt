# CmakeLists.txt in DemoUI

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Qml)
set(QT_USE_QTSVG TRUE)
set(QT_USE_QTSCRIPT TRUE)

include_directories(${Cornucopia_SOURCE_DIR}/Cornucopia)
include_directories(${Cornucopia_SOURCE_DIR}/DemoUI)

file(GLOB DemoUI_CPP "*.cpp")
file(GLOB DemoUI_H "*.h")
file(GLOB DemoUI_UI "*.ui")

qt_wrap_ui(DemoUI_H ${DemoUI_UI})

list(APPEND DemoUI_Sources ${DemoUI_CPP} ${DemoUI_H})

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(DemoUI ${DemoUI_Sources})

# This is for moc files
set_target_properties(DemoUI PROPERTIES AUTOMOC TRUE)

target_link_libraries(DemoUI PRIVATE Cornucopia Qt6::Core Qt6::Gui Qt6::Widgets
                                     Qt6::Qml)

install(TARGETS DemoUI RUNTIME DESTINATION bin)

if(WIN32)
  install(
    FILES "${QT_BINARY_DIR}/QtCore${QT_VERSION_MAJOR}.dll"
          "${QT_BINARY_DIR}/QtGui${QT_VERSION_MAJOR}.dll"
          "${QT_BINARY_DIR}/QtScript${QT_VERSION_MAJOR}.dll" DESTINATION bin)
endif()
