#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>
#include <mfa/mfa.hpp>
#include "opts.h"
#include "block.hpp"


namespace mfa_extend
{

    // template<typename T>
    // void recover_mfa_selected(
    //         const mfa::MFA<T>* mfa,
    //         const Block<T>* b,               // mfa data model
    //         VectorXi&            span,       //span index should be in [p, n]
    //         VectorX<T>&   point,                  // parameters of point to decode
    //         VectorX<T>&             weights, //weights for control points
    //         VectorX<T>&         cpt,const VectorXi&     derivs = VectorXi())
    // {
    //     VectorX<T> domain_range = b->core_maxs - b->core_mins;
    //     VectorX<T> param = (point-b->core_mins).cwiseQuotient(domain_range);

    //     mfa->DecodePt_selectd_span(*mfa_data,span, param,derivs,weights,cpt);
    //     for(int i=0;i<derivs.size();++i)
    //     {
    //         if(derivs[i]>0)
    //         {
    //             cpt = cpt/pow(domain_range[i],derivs[i]);
    //         }
    //     }
    // }
    
    // template<typename T>
    // void recover_mfa_selected(
    //         mfa::MFA<T>* mfa,
    //         mfa::MFA_Data<T>*  mfa_data,               // mfa data model
    //         VectorXi&            span,       //span index should be in [p, n]
    //         VectorX<T>&   point,                  // parameters of point to decode
    //                                                     // pass size-0 vector if unused
    //         //MatrixX<T>&             ctrl_pts,   //control points of first derivative
    //         VectorX<T>&             weights, //weights for control points
    //         VectorX<T>&         cpt,
    //         const VectorX<T>&         domain_min,
    //         const VectorX<T>&         domain_range)
    // {
    //     VectorX<T> param = (point-domain_min).cwiseQuotient(domain_range);

    //     mfa->DecodePt_selectd_span(*mfa_data,span, param,weights,cpt);
    // }

    template<typename T>
    void recover_mfa(
            const Block<T>* b,               // mfa data model
            const VectorX<T>&   point,                  // parameters of point to decode
            VectorX<T>&         cpt,const VectorXi&     derivs = VectorXi())
    {
        VectorX<T> domain_range = b->core_maxs - b->core_mins;
        VectorX<T> param = (point-b->core_mins).cwiseQuotient(domain_range);
        b->mfa->DecodeVar(0,param,cpt,derivs);

        if (derivs.size() > 0)
        {
            for(int i=0;i<derivs.size();++i)
            {
                if(derivs[i]>0)
                {
                    cpt = cpt/pow(domain_range[i],derivs[i]);
                }
            }
        }
    }

}