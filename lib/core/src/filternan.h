#ifndef FILTERNAN_H
#define FILTERNAN_H

namespace IPP
{

template<typename T>
class FilterNaN
{
public:
    FilterNaN()
        : fallBackVal(0.0)
    {

    }

    FilterNaN(const T& fallBackVal)
        : fallBackVal(fallBackVal)
    {

    }

    T operator()(T& input)
    {
        if(input == input)
        {
            return input;
        }
        else
        {
            std::cout << input << "\t" << fallBackVal << std::endl;
            return fallBackVal;
        }
    }

private:
    const T fallBackVal;
};


}


#endif // FILTERNAN_H
