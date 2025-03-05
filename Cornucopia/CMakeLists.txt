# Cornucopia library CMakeLists.txt

file(GLOB Cornucopia_CPP "*.cpp")
file(GLOB Cornucopia_H "*.h")
# list(APPEND Cornucopia_Sources ${Cornucopia_CPP} ${Cornucopia_H})

add_library(Cornucopia STATIC ${Cornucopia_CPP})

install(TARGETS Cornucopia ARCHIVE DESTINATION lib)

# install( DIRECTORY . DESTINATION include FILES_MATCHING PATTERN "*.h")

# Install headers
install(FILES ${Cornucopia_HEADERS} DESTINATION include)
