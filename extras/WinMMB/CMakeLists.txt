cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(WinMMB)

add_definitions("-DWIN32_LEAN_AND_MEAN")

if (MSVC)
    add_definitions("-D_MSVC_COMPILER")
else()
    add_definitions("-D_UNICODE -DUNICODE")
endif()

if (MINGW)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -pedantic -Wall -Wextra")
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/src")

set(WinMMB_SRCS
    "src/main.c"
    "src/Resource.rc"
)

add_executable(WinMMB WIN32 ${WinMMB_SRCS})
target_link_libraries(WinMMB comctl32)
