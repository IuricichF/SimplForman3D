cmake_minimum_required(VERSION 2.8)

file(GLOB SOURCE_LIBFORMAN "*.h" "*.cpp")
add_library(libFORMAN STATIC ${SOURCE_LIBFORMAN})

TARGET_LINK_LIBRARIES(libFORMAN libIASTAR)
