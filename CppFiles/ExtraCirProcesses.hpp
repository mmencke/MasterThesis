//
//  ExtraCirProcesses.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 25/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef ExtraCirProcesses_hpp
#define ExtraCirProcesses_hpp

//#include <iostream>

#include <ql/stochasticprocess.hpp>
#include <ql/processes/eulerdiscretization.hpp>

namespace QuantLib {
    
    //! CoxIngersollRoss process class
    /*! This class describes the CoxIngersollRoss process governed by
     \f[
     dx = a (r - x_t) dt + \sqrt{x_t}\sigma dW_t.
     \f]
     \ingroup processes
     */
    class CirProcess : public StochasticProcess1D {
    public:
        enum Discretization { None,
            FullTruncation,
            QuadraticExponential,
            Exact
        };
        CirProcess(Real speed,
                                Volatility vol,
                                Real x0 = 0.0,
                                Real level = 0.0,
                                Discretization d = None);
        //@{
        Real drift(Time t, Real x) const override;
        Real diffusion(Time t, Real x) const override;
        Real expectation(Time t0, Real x0, Time dt) const override;
        Real stdDeviation(Time t0, Real x0, Time dt) const override;
        //@}
        Real x0() const override;
        Real speed() const;
        Real volatility() const;
        Real level() const;
        Real variance(Time t0, Real x0, Time dt) const override;
        Real evolve (Time t0,
                     Real x0,
                     Time dt,
                     Real dw) const override;
    private:
        Real x0_, speed_, level_;
        Volatility volatility_;
        Discretization discretization_;
    };
    
    
    
    class CirHelperProcess : public StochasticProcess1D {
        public:
            CirHelperProcess( Real k, Real sigma, Real y0, Real theta)
            : y0_(y0), theta_(theta), k_(k), sigma_(sigma) {
                discretization_ =
                ext::shared_ptr<discretization>(new EulerDiscretization);
            }
        
            Real x0() const override { return y0_; }
            Real drift(Time, Real y) const override {
                Real drift_adj;
                if(std::fabs(y)<0.001) {
                    drift_adj=0; //make sure that the process does not explode when y approaches 0
                } else {
                    drift_adj=(0.5*theta_*k_ - 0.125*sigma_*sigma_)/y
                    - 0.5*k_*y;
                }
                return drift_adj;
            }
            Real diffusion(Time, Real) const override { return 0.5 * sigma_; }
        
        private:
            Real y0_, theta_, k_, sigma_;
    };

    
}

#endif /* ExtraCirProcesses_hpp */
