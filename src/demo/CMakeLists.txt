set(target softpositdemo)
add_executable(${target} demo.cpp)
target_link_libraries(${target} softposit)
