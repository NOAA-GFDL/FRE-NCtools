# https://approvaltestscpp.readthedocs.io/en/latest/generated_docs/CMakeIntegration.html#make-cmake-clone-approvaltests-cpp-and-catch2

# Needs CMake 3.14 or above
Include(FetchContent)

# Names of targets added here:
#   ghc_filesystem Catch2::Catch2 ApprovalTests::ApprovalTests fmt::fmt
#TODO: update git tags to latest stables?

# -------------------------------------------------------------------
FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0)
FetchContent_MakeAvailable(Catch2)
###        v2.13.4? devel? 3.3.3 (latest of 3.)

# -------------------------------------------------------------------
FetchContent_Declare(ApprovalTests
        GIT_REPOSITORY https://github.com/approvals/ApprovalTests.cpp.git
        GIT_TAG master)
FetchContent_MakeAvailable(ApprovalTests)

# -------------------------------------------------------------------
## The format libray was instead  brought in with the rest of the C++23 features,
## and added to the std namespace. See the toplevel CMakeLists file.
#FetchContent_Declare(fmt
#        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
#        GIT_TAG  7.1.3)
#FetchContent_MakeAvailable(fmt)

#FetchContent_Declare(fmt
#        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
#       GIT_TAG 10.0.0
#        SYSTEM
#        OVERRIDE_FIND_PACKAGE)
#FetchContent_MakeAvailable(fmt)

# -------------------------------------------------------------------
FetchContent_Declare(filesystem
        GIT_REPOSITORY https://github.com/gulrak/filesystem.git
        GIT_TAG v1.5.0)
FetchContent_MakeAvailable(filesystem)
