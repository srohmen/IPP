#ifndef PLBSERIALIZATION_H
#define PLBSERIALIZATION_H


namespace boost
{
namespace serialization
{


template<class Archive>
static void serialize(Archive & ar, plb::Dot2D& dot, const unsigned int version)
{
    ar & dot.x;
    ar & dot.y;
}

template<class Archive>
static void serialize(Archive & ar, plb::DotList2D& list, const unsigned int version)
{
    ar & list.dots;
}

template<class Archive>
static void serialize(Archive & ar, plb::Dot3D& dot, const unsigned int version)
{
    ar & dot.x;
    ar & dot.y;
    ar & dot.z;
}

template<class Archive>
static void serialize(Archive & ar, plb::DotList3D& list, const unsigned int version)
{
    ar & list.dots;
}

}
}


#endif // PLBSERIALIZATION_H
