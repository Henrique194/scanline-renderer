set(SOURCE_FILES
    input.c
    input.h
    main.c
    palette.c
    palette.h
    video.c
    video.h
)
set(LIBS ${SDL2_LIBRARIES})

if(WIN32)
    add_executable(${PACKAGE_TARNAME} WIN32 ${SOURCE_FILES})
else()
    add_executable(${PACKAGE_TARNAME} ${SOURCE_FILES})
endif()

target_include_directories(${PACKAGE_TARNAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PACKAGE_TARNAME} ${LIBS})

if(MSVC)
    set_target_properties(${PACKAGE_TARNAME} PROPERTIES LINK_FLAGS "/MANIFEST:NO")
endif()
