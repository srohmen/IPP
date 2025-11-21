#ifndef PALABOSFIELDDECOMPOSITION_H
#define PALABOSFIELDDECOMPOSITION_H

#include "plbtypededuction.h"
#include "fielddecomposition.h"
#include "palabosmpiutils.h"
#include "dimsizecalc.h"

namespace IPP
{

namespace ConvertHelper
{

typedef std::vector<plb::plint> PlbIntVec;


template<size_t dim>
struct ConvertToBox;


template<>
struct ConvertToBox<2>
{
    static IPPBox3DLong convert(PlbIntVec::const_iterator begin)
    {
        return IPPBox3DLong( {{ *begin, *(begin+2), 0 }}, {{ *(begin+1), *(begin+3), 0 }} );
    }

    static IPPBox3DLong convert(const plb::Box2D &box)
    {
        return IPPBox3DLong( {{ box.x0, box.y0, 0 }}, {{ box.x1, box.y1, 0 }} );
    }
};


template<>
struct ConvertToBox<3>
{
    static IPPBox3DLong convert(PlbIntVec::const_iterator begin)
    {
        return IPPBox3DLong(
        {{ *begin, *(begin+2), *(begin+4) }},
        {{ *(begin+1), *(begin+3), *(begin+5) }}
                    );
    }

    static IPPBox3DLong convert(const plb::Box3D &box)
    {
        return IPPBox3DLong(
        {{ box.x0, box.y0, box.z0 }},
        {{ box.x1, box.y1, box.z1 }}
                    );
    }
};


template<size_t dim>
void convertToBoxes(const std::vector<plb::plint>& allGlobalIndices,
                    std::vector<IPPBox3DLong>& localDomains)
{
    constexpr size_t nComp = dim * 2;
    for(PlbIntVec::const_iterator it = allGlobalIndices.begin();
        it < allGlobalIndices.end(); it += nComp)
    {
        localDomains.push_back(ConvertToBox<dim>::convert(it));
    }
}

inline plb::Box2D operator+(const plb::Box2D& box, const plb::Dot2D& dot)
{
    return plb::Box2D(box.x0 + dot.x,
                      box.x1 + dot.x,
                      box.y0 + dot.y,
                      box.y1 + dot.y);
}

inline plb::Box3D operator+(const plb::Box3D& box, const plb::Dot3D& dot)
{
    return plb::Box3D(box.x0 + dot.x,
                      box.x1 + dot.x,
                      box.y0 + dot.y,
                      box.y1 + dot.y,
                      box.z0 + dot.z,
                      box.z1 + dot.z);
}

}

template<typename T, typename Field>
class RetrievePalabosDecomposedDomains :
        public PlbTypeDeduction::GetBoxProcessingFunc<T, Field>::value
{
public:
    using BaseClass = typename PlbTypeDeduction::GetBoxProcessingFunc<T, Field>::value;
    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;


    RetrievePalabosDecomposedDomains(size_t& localNumberCells,
                                     IPPBox3DLong& ownDomain,
                                     std::vector<IPPBox3DLong> & localDomains,
                                     std::vector<int>& displs,
                                     std::vector<int>& sndCounts)
        : localNumberCells(localNumberCells),
          ownDomain(ownDomain),
          localDomains(localDomains),
          displs(displs),
          sndCounts(sndCounts)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
    }

    virtual RetrievePalabosDecomposedDomains<T, Field>* clone() const
    {
        return new RetrievePalabosDecomposedDomains<T, Field>(*this);
    }


    virtual void process(Box domain, InputField& decomposedMultiData)
    {
        using namespace ConvertHelper;
        constexpr size_t dim = PlbTypeDeduction::GetFieldDim<T, Field>::value;

        const auto offset = decomposedMultiData.getLocation();
        localNumberCells = domain.nCells();
        ownDomain = ConvertToBox<dim>::convert(domain + offset);

        const bool isMainProc = plb::global::mpi().isMainProcessor();
        const size_t nProcesses = plb::global::mpi().getSize();

        const PalabosMPIUtils::ProcInfo procInfo(isMainProc, nProcesses);

        std::vector<plb::plint> allGlobalIndices;
        PalabosMPIUtils::allGatherAllIndices(procInfo, offset,
                                          domain, allGlobalIndices);

        if(isMainProc)
        {
            using SizeCalc = DimSizeCalc<dim>;
            PalabosMPIUtils::transferStrides<SizeCalc>(allGlobalIndices, nProcesses,
                                                       sndCounts, displs);


        }

        convertToBoxes<dim>(allGlobalIndices, localDomains);
    }


private:
    size_t& localNumberCells;
    IPPBox3DLong& ownDomain;
    std::vector<IPPBox3DLong>& localDomains;
    std::vector<int>& displs;
    std::vector<int>& sndCounts;

};


template<typename T, typename MultiField>
class PalabosFieldDecomposition : public FieldDecomposition
{
public:
    PalabosFieldDecomposition(MultiField& field)
    {
        MPI_CHECK_SYNC;

        const auto box = field.getBoundingBox();
        globalSize = createSizeVec(box);
        using DomainRecieve = RetrievePalabosDecomposedDomains<T, MultiField>;
        plb::applyProcessingFunctional(new DomainRecieve(localNumberCells,
                                                         ownDomain,
                                                         domains,
                                                         displs,
                                                         sndCounts),
                                       field.getBoundingBox(), field);

        MPI_CHECK_SYNC;
    }

private:
    static IPPVector3DLong createSizeVec(const plb::Box2D& box)
    {
        return IPPVector3DLong({{ (box.x1 - box.x0) + 1,
                                  (box.y1 - box.y0) + 1,
                                  1
                                }});
    }

    static IPPVector3DLong createSizeVec(const plb::Box3D& box)
    {
        return IPPVector3DLong({{ (box.x1 - box.x0) + 1,
                                  (box.y1 - box.y0) + 1,
                                  (box.z1 - box.z0) + 1
                                }});

    }
};

}

#endif // PALABOSFIELDDECOMPOSITION_H
