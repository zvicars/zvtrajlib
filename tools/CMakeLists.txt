add_library (TOOLS SHARED
Stacktrace.hpp
Stacktrace.cpp
Assert.hpp
GenericFactory.hpp
InputParser.cpp
InputParser.hpp
StringTools.cpp
StringTools.hpp
)
target_include_directories(TOOLS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})