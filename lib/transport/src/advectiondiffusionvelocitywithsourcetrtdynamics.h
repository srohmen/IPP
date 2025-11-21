#ifndef ADVECTIONDIFFUSIONVELOCITYWITHSOURCETRTDYNAMICS_H
#define ADVECTIONDIFFUSIONVELOCITYWITHSOURCETRTDYNAMICS_H

#include <palabos/complexDynamics/trtDynamics.h>

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionVelocityWithSourceTRTDynamics : public plb::TRTdynamics<T, Descriptor>
{
public:
    AdvectionDiffusionVelocityWithSourceTRTDynamics(const T& omega);

    virtual AdvectionDiffusionVelocityWithSourceTRTDynamics<T,Descriptor>* clone() const;

    virtual int getId() const;

    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics)
    {
        plb::Array<T,Descriptor<T>::q> eq;
        // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
        // So we allocate the index-zero memory location, and waste some memory
        // for convenience.
        plb::Array<T,Descriptor<T>::q/2+1> eq_plus, eq_minus, f_plus, f_minus;


        plb::Array<T,Descriptor<T>::d> j;
        T rhoBar;
        plb::momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
        rhoBar /= poros;
        j[0] = 0;
        j[1] = 0;
        T jSqr = normSqr(j);
        T invRho = Descriptor<T>::invRho(rhoBar);

        equilibrium<T, Descriptor<T>>(rhoBar, j, jSqr, eq);



        for (plb::plint i=1; i<=Descriptor<T>::q/2; ++i)
        {
            const plb::plint iOpp = i + Descriptor<T>::q / 2;
            eq_plus[i]  = 0.5*(eq[i] + eq[iOpp]);
            eq_minus[i] = 0.5*(eq[i] - eq[iOpp]);
            f_plus[i]   = 0.5*(cell[i] + cell[iOpp]);
            f_minus[i]  = 0.5*(cell[i] - cell[iOpp]);
        }



        // tau_s = 0.5 + (MagicPara/(tau_a(i, j)-0.5))
        // omega_a = 1.d0 / tau_a (i, j)
        // omega_s = 1.d0 / tau_s
        //    const T w2 = cPhi * w;
        //    const T w1 = 1.0 - 2.0 * cPhi;

        const T sMinus = this->getOmega();
        const T sPlus = s_magic / (1.0/sMinus - 0.5) + 0.5;

        // cell[0] += -sPlus * cell[0] + sPlus * eq[0];
        cell[0] += sPlus * (eq[0] - cell[0]);

        for (plb::plint i=1; i<=Descriptor<T>::q/2; ++i)
        {
            cell[i] += -sPlus * (f_plus[i] - eq_plus[i]) - sMinus * (f_minus[i] - eq_minus[i]);
            cell[i+Descriptor<T>::q/2] += -sPlus * (f_plus[i] - eq_plus[i]) + sMinus * (f_minus[i] - eq_minus[i]);
        }


        if (cell.takesStatistics())
        {
            gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho );
        }
    }

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(plb::Cell<T,Descriptor>& cell, T rhoBar,
                                 plb::Array<T,Descriptor<T>::d> const& j,
                                 T thetaBar, plb::BlockStatistics& stat)
    {
        PLB_ASSERT(false);
    }

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plb::plint iPop, T rhoBar, plb::Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const
    {
        if(iPop != 0)
        {
            const T f = equilibrium<T, Descriptor<T>>(iPop, rhoBar, j, jSqr);
            return f;
        }
        else
        {
            plb::Array<T, Descriptor<T>::q> eq;
            equilibrium<T, Descriptor<T>>(rhoBar, j, jSqr, eq);
            return eq[0];
        }
    }

    virtual T computeDensity(plb::Cell<T,Descriptor> const& cell) const
    {
        T rhoBar = plb::momentTemplates<T,Descriptor>::get_rhoBar(cell) / poros;
        return rhoBar;
    }


private:
    static const T s_magic;
    static int id;
};

}

#endif // ADVECTIONDIFFUSIONVELOCITYWITHSOURCETRTDYNAMICS_H
