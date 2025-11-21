#ifndef PLB_HASH_H
#define PLB_HASH_H

namespace std
{

template <>
struct hash<plb::Dot2D>
{
    typedef plb::Dot2D argument_type;
    typedef std::size_t result_type;

    result_type operator()(const plb::Dot2D & t) const
    {
        std::size_t val { 0 };
        boost::hash_combine(val, t.x);
        boost::hash_combine(val, t.y);
        return val;
    }
};

template <>
struct hash<plb::Dot3D>
{
    typedef plb::Dot3D argument_type;
    typedef std::size_t result_type;

    result_type operator()(const plb::Dot3D & t) const
    {
        std::size_t val { 0 };
        boost::hash_combine(val, t.x);
        boost::hash_combine(val, t.y);
        boost::hash_combine(val, t.z);
        return val;
    }
};

}

#endif // PLB_HASH_H
