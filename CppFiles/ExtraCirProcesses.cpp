//
//  ExtraCirProcesses.cpp
//  Master Thesis
//
//  Created by Magnus Mencke on 25/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#include "ExtraCirProcesses.hpp"


#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <ql/math/distributions/normaldistribution.hpp>

namespace QuantLib {
    Real CirProcess::x0() const {
        return x0_;
    }
    
    Real CirProcess::speed() const {
        return speed_;
    }
    
    Real CirProcess::volatility() const {
        return volatility_;
    }
    
    Real CirProcess::level() const {
        return level_;
    }
    
    Real CirProcess::drift(Time, Real x) const {
        return speed_ * (level_ - x);
    }
    
    Real CirProcess::diffusion(Time, Real) const {
        return volatility_;
    }
    
    Real CirProcess::expectation(Time, Real x0,
                                        Time dt) const {
        return level_ + (x0 - level_) * std::exp(-speed_*dt);
    }
    
    Real CirProcess::stdDeviation(Time t, Real x0,
                                         Time dt) const {
        return std::sqrt(variance(t,x0,dt));
    }

    //Full truncation scheme (as in Brigo, Morini and Pallavicini)
    //Inspired by Heston implementation
    //  see Lord, R., R. Koekkoek and D. van Dijk (2006),
    // "A Comparison of biased simulation schemes for
    //  stochastic volatility models",
    // Working Paper, Tinbergen Institute
    Real CirProcess::evolve (Time t0,
                                    Real x0,
                                    Time dt,
                                    Real dw) const {
        Real resultTrunc;
        switch (discretization_) {
            case None: {
                resultTrunc=apply(expectation(t0,x0,dt),stdDeviation(t0,x0,dt)*dw);
                break;
            }
            case FullTruncation: {
                Real x0_trunc = x0>0.0 ? x0 : 0.0;
                
                Real result = apply( expectation(t0, x0_trunc, dt),stdDeviation(t0,x0_trunc,dt)*dw);
                
                resultTrunc = result>0.0 ? result : 0.0;
                
                break;
            }
            case QuadraticExponential: {
                // for details of the quadratic exponential discretization scheme
                // see Leif Andersen,
                // Efficient Simulation of the Heston Stochastic Volatility Model
                const Real ex = std::exp(-speed_*dt);
                
                const Real m  =  level_+(x0-level_)*ex;
                const Real s2 =  x0*volatility_*volatility_*ex/speed_*(1-ex)
                + level_*volatility_*volatility_/(2*speed_)*(1-ex)*(1-ex);
                const Real psi = s2/(m*m);
                
                if (psi <= 1.5) {
                    const Real b2 = 2/psi-1+std::sqrt(2/psi*(2/psi-1));
                    const Real b  = std::sqrt(b2);
                    const Real a  = m/(1+b2);
                    
                    resultTrunc = a*(b+dw)*(b+dw);
                }
                else {
                    const Real p = (psi-1)/(psi+1);
                    const Real beta = (1-p)/m;
                    
                    const Real u = CumulativeNormalDistribution()(dw);
                    
                    resultTrunc = ((u <= p) ? 0.0 : std::log((1-p)/(1-u))/beta);
                }
                break;
            }
            case Exact: {
                CumulativeNormalDistribution dwDist; //despite the name, dw is standard normal
                Real uniform= dwDist(dw); //transforming normal to uniform
                
                Real c=(4*speed_)/(volatility_*volatility_*(1-std::exp(-speed_*dt)));
                Real nu=(4*speed_*level_)/(volatility_*volatility_);
                Real eta=c*x0*std::exp(-speed_*dt);
                
                //we had some strange errors using QuantLib's Chi2 distribution, so we use the one from boost
                boost::math::non_central_chi_squared_distribution<double,  boost::math::policies::policy<>> chi2(nu, eta);
                //InverseNonCentralCumulativeChiSquareDistribution chi2(nu, eta);
                
                resultTrunc = quantile(chi2,uniform)/c;
            }
        }
        
        return resultTrunc;
    }
    
    
    CirProcess::CirProcess(Real speed,
                                                     Volatility vol,
                                                     Real x0,
                                                     Real level,
                                                     Discretization d)
    : x0_(x0), speed_(speed), level_(level), volatility_(vol), discretization_(d) {
        QL_REQUIRE(volatility_ >= 0.0, "negative volatility given");
    }
    
    Real CirProcess::variance(Time, Real, Time dt) const {
        Real exponent1 = std::exp(-speed_ * dt);
        Real exponent2 = std::exp(-2 * speed_ * dt);
        Real fraction = (volatility_ * volatility_) / speed_;
        
        return x0_ * fraction * (exponent1 - exponent2) + level_ * fraction * (1 - exponent1) * (1 - exponent1);
    }
    
}
