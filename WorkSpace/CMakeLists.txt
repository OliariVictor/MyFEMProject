file(GLOB headers *.h Bases/*.h)
file(GLOB sources *.cpp Bases/*.cpp)

include_directories(./Bases)

message("sources: ${sources}")

set(Mymain main.cpp ${headers} ${sources})

add_executable(HdivEx ${Mymain})

target_link_libraries(HdivEx pz)
