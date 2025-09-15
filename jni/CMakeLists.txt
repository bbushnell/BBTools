#Project settings
cmake_minimum_required (VERSION 2.8) 
project(BBToolsJNI) 

#Find required packages
find_package(Java REQUIRED)
find_package(JNI REQUIRED)
if(NOT EXISTS ${Java_JAVAH_EXECUTABLE})
    FIND_PROGRAM(Java_JAVAH_EXECUTABLE
        NAMES javah
        HINTS ${_JAVA_HINTS}
        PATHS ${_JAVA_PATHS}
    )
    if(EXISTS ${Java_JAVAH_EXECUTABLE})
        message(STATUS "Found Java_JAVAH_EXECUTABLE=${Java_JAVAH_EXECUTABLE}")
    else()
        message(FATAL_ERROR "Java generator 'javah' could not be found. Ensure it is installed and in your PATH.")
    endif()
endif()

#Set C_FLAGS on Unix
if(UNIX OR APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
    set(CMAKE_C_FLAGS "-O3 -std=c99")
endif()

#Set the filename prefixes
set(NAME1 "BandedAlignerJNI")
set(NAME2 "MultiStateAligner11tsJNI")
set(NAME3 "BBMergeOverlapper")
set(NAME4 "IceCreamAlignerJNI")

#Set names of source and header files
set(SOURCES "${NAME1}.c" "${NAME2}.c" "${NAME3}.c" "${NAME4}.c")
set(HEADERS "align2_${NAME1}.h" "align2_${NAME2}.h" "jgi_${NAME3}.h" "jgi_${NAME4}.h")

#Generate header files
set(_stubDir "${CMAKE_CURRENT_BINARY_DIR}")
set(_classDir "${CMAKE_BINARY_DIR}")
add_custom_command(
    OUTPUT align2_${NAME1}.h
    COMMAND ${Java_JAVAH_EXECUTABLE} -verbose -classpath ${_classDir} -d ${_stubDir} -jni align2.${NAME1}
    )
add_custom_command(
    OUTPUT align2_${NAME2}.h
    COMMAND ${Java_JAVAH_EXECUTABLE} -verbose -classpath ${_classDir} -d ${_stubDir} -jni align2.${NAME2}
    )
add_custom_command(
    OUTPUT jgi_${NAME3}.h
    COMMAND ${Java_JAVAH_EXECUTABLE} -verbose -classpath ${_classDir} -d ${_stubDir} -jni jgi.${NAME3}
    )
add_custom_command(
    OUTPUT jgi_${NAME3}.h
    COMMAND ${Java_JAVAH_EXECUTABLE} -verbose -classpath ${_classDir} -d ${_stubDir} -jni jgi.${NAME4}
    )

#Generate bbtoolsjni.jnilib
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${JNI_INCLUDE_DIRS}) 
add_library(bbtoolsjni SHARED ${SOURCES} ${HEADERS})

