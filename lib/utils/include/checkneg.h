#ifndef CHECKNEG_H
#define CHECKNEG_H


namespace IPP
{

template<typename T>
class CheckNeg
{
public:
    CheckNeg(bool& hasNeg)
        : m_hasNeg(hasNeg)
    {

    }

    T operator()(const T& input)
    {
        if(input < 0)
        {
            m_hasNeg = true;
            assert(false);
        }

        return input;
    }

private:
    bool& m_hasNeg;
};


}


#endif // CHECKNEG_H
