cmake_minimum_required(VERSION 3.15)

project(formula_weird_engine
    VERSION 0.1.0
    DESCRIPTION "C++17 explicit-dynamics physics engine"
    LANGUAGES CXX
)

option(FORMULA_WEIRD_ENABLE_OPENMP "Enable OpenMP where supported" ON)
option(FORMULA_WEIRD_STRICT_WARNINGS "Enable strict compiler warnings" ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if(FORMULA_WEIRD_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

set(FORMULA_WEIRD_SOURCES
    src/constitutive/constitutive_model.cpp
    src/thermo/thermodynamics.cpp
    src/thermo/thermodynamic_consistency.cpp
    src/failure/failure.cpp
    src/rigid/rigid_body_sync.cpp
    src/debris/debris.cpp
    src/stability/stability.cpp
    src/time/time_integration.cpp
    src/contact/contact_mapping.cpp
)

add_executable(formula_weird_engine
    src/main.cpp
    ${FORMULA_WEIRD_SOURCES}
)

target_compile_features(formula_weird_engine PRIVATE cxx_std_17)

target_include_directories(formula_weird_engine
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(formula_weird_engine
    PRIVATE
        Eigen3::Eigen
)

if(FORMULA_WEIRD_ENABLE_OPENMP)
    target_link_libraries(formula_weird_engine PRIVATE OpenMP::OpenMP_CXX)
endif()

set_target_properties(formula_weird_engine PROPERTIES
    CXX_EXTENSIONS OFF
)

if(FORMULA_WEIRD_STRICT_WARNINGS)
    if(MSVC)
        target_compile_options(formula_weird_engine PRIVATE /W4 /permissive-)
    else()
        target_compile_options(formula_weird_engine PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow
        )
    endif()
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(remaining_modules_tests
        tests/remaining_modules_tests.cpp
        src/debris/debris.cpp
        src/stability/stability.cpp
        src/time/time_integration.cpp
        src/contact/contact_mapping.cpp
    )
    target_compile_features(remaining_modules_tests PRIVATE cxx_std_17)
    target_include_directories(remaining_modules_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(remaining_modules_tests PRIVATE Eigen3::Eigen)
    if(FORMULA_WEIRD_STRICT_WARNINGS)
        if(MSVC)
            target_compile_options(remaining_modules_tests PRIVATE /W4 /permissive-)
        else()
            target_compile_options(remaining_modules_tests PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
        endif()
    endif()
    add_test(NAME remaining_modules_tests COMMAND remaining_modules_tests)
endif()
