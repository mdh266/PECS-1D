# CMAKE FILE 
##########################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# set the ./run directory to be the directory where to put 
# the executables
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_INSTALL_PREFIX_PATH
		./run/)

#Executable Name
PROJECT(PROJECT)


#find the eigen directory
FIND_PATH(EIGEN_INCLUDE_DIR eigen3/signature_of_eigen3_matrix_library)
SET(EIGEN_DIR ${EIGEN_INCLUDE_DIR}/eigen3/)

#find GSL
INCLUDE(FindGSL)
FIND_PACKAGE(GSL REQUIRED)

#find BOOST
FIND_PACKAGE(Boost REQUIRED)


#find OpenMP (optional)
FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()


SET(HEADERS "include")

# set the include directories
SET(INCLUDE_DIRS
	${Boost_INCLUDE_DIRS}
	${GSL_INCLUDE_DIR}
	${EIGEN_DIR}
	${HEADERS}
	)

# set the libraries to link to
SET(LIBS 
	 ${Boost_LIBRARIES}
	 ${GSL_LIBRARIES}
	 )

# set the project_source directories
SET(PROJECT_SRCS
	 ./src/quadrule.cpp
	 ./src/readAndPrint.cpp
	 ./ConcentrationProfiles/ConcentrationProfiles.cpp
	 ./src/Utilities.cpp
	 ./src/ddpComputeDeltaT.cpp
	 ./src/problemSetUp.cpp
	 ./src/basis.cpp
	 ./src/VandeMondeSetUp.cpp
	 ./src/matrixSetUp.cpp
	 ./src/vectorSetUp.cpp
	 ./src/makeProperties.cpp
	 ./src/testUtilities.cpp
	./src/carrier.cpp
	 ./src/poisson.cpp
	 ./src/dopingprofile.cpp
	 ./src/test.cpp
	 ./src/ddpPrintState.cpp
	)


SET(SOLARCELL_SRCS
	 ./src/ddpTimeStepping.cpp
	 ./src/ddpRecombination.cpp
	)


SET(TEST_SRCS
	./src/test.cpp
	)

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

ADD_EXECUTABLE(solar_cell_app 
							./src/CTU_IMEX.cpp 
							${PROJECT_SRCS} 
							${SOLARCELL_SRCS}
							)
ADD_EXECUTABLE(mesfet 
							./src/Implicit_Mesfet.cpp 
							${PROJECT_SRCS}
							)

ADD_EXECUTABLE(test_System
							./src/test_System.cpp
							${PROJECT_SRCS}
							${TEST_SRCS}
							)

ADD_EXECUTABLE(L2_Schwarz 
							./src/L2_Schwarz_Errors.cpp
							${PROJECT_SRCS}
							${SOLARCELL_SRCS}
							)

TARGET_LINK_LIBRARIES(solar_cell_app ${LIBS})
TARGET_LINK_LIBRARIES(mesfet ${LIBS})
TARGET_LINK_LIBRARIES(test_System ${LIBS})
TARGET_LINK_LIBRARIES(L2_Schwarz ${LIBS})


