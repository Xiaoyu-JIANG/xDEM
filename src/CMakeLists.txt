add_library(
    xDEM STATIC
    BaseContactModel.h
    BaseParticle.cpp
    BaseParticle.h
    BaseWorld.cpp
    BaseWorld.h
    CircularParticle.cpp
    CircularParticle.h
    ClusterLabelling.cpp
    ClusterLabelling.h
    ContactGeometry.h
    ContactInformation.h
    ContactModelList.h
    ContactModel_LinearDissipativeFrictional.cpp
    ContactModel_LinearDissipativeFrictional.h
    ContactModel_UniversalGravity.cpp
    ContactModel_UniversalGravity.h
    EllipticParticle.cpp
    EllipticParticle.h
    GluedParticle.cpp
    GluedParticle.h
    LeesEdwardsBoundaryWorld.cpp
    LeesEdwardsBoundaryWorld.h
    LongRangeForceWorld.cpp
    LongRangeForceWorld.h
    MathMatrix2d.h
    MathVector.h
    PeriodicBoundaryWorld.cpp
    PeriodicBoundaryWorld.h
    PeriodicControlWorld.cpp
    PeriodicControlWorld.h
    PointParticle.cpp
    PointParticle.h
    QuadronAnalysis.cpp
    QuadronAnalysis.h
    SegmentObject.cpp
    SegmentObject.h
)


find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()