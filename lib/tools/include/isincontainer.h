#ifndef ISINCONTAINER_H
#define ISINCONTAINER_H

template<typename ContainerType>
class IsInContainer
{
public:
    IsInContainer(const ContainerType& container)
        : container(container)
    {

    }

    bool operator()(const typename ContainerType::value_type& elem)
    {
        return container.find(elem) != container.end();
    }


private:
    const ContainerType& container;
};

#endif // ISINCONTAINER_H
