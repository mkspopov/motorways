cmake_minimum_required(VERSION 3.21)
project(motorways)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")

add_executable(motorways main.cpp)

find_package(SFML REQUIRED COMPONENTS graphics audio window system)

target_link_libraries(
    ${PROJECT_NAME} PUBLIC
    sfml-graphics
    sfml-audio
    sfml-window
    sfml-system
    sfml-network
    pthread
)
