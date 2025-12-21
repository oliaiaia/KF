cmake_minimum_required(VERSION 3.10)
project (KalmanFilter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(
  include
  ${CMAKE_SOURCE_DIR}/include
)

add_subdirectory(matplotplusplus)

find_package(Eigen3 3.3 REQUIRED )
add_executable(kalmanFilter src/exmpl0.cpp src/KF.cpp)
target_link_libraries(kalmanFilter Eigen3::Eigen matplot)

add_executable(EKF2sources src/exmpl2.cpp src/EKF.cpp)
target_link_libraries(EKF2sources Eigen3::Eigen matplot)

# add_executable(EKF src/exmpl1.cpp src/EKF_11.cpp)
# target_link_libraries(EKF Eigen3::Eigen matplot)
