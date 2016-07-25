/*
 * Copyright (c) 2015, Jonathan Ventura
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef POLYNOMIAL_INTERNAL_HPP
#define POLYNOMIAL_INTERNAL_HPP

#include <Eigen/Core>
#include <vector>
#include <queue>

namespace PolynomialVentura
{
    namespace Internal
    {
        
        template<int a,int b> struct max{ static const int value=(a<b)?b:a; };
        template<int a,int b> struct min{ static const int value=(a>b)?b:a; };
        
        template <int length>
        inline
        Eigen::Map< Eigen::Matrix<double,1,length> > vecmap(double *data)
        {
            return Eigen::Map< Eigen::Matrix<double,1,length> >(data);
        }
        
        template <int length>
        inline
        Eigen::Map< const Eigen::Matrix<double,1,length> > vecmap(const double *data)
        {
            return Eigen::Map< const Eigen::Matrix<double,1,length> >(data);
        }
        
        template <int deg1,int deg2,int k>
        struct
        _PolyConv
        {
            static void
            compute( double *coefout, const double *coef1, const double *coef2 )
            {
                coefout[k] =
                (
                 vecmap<min<k,deg1>::value-max<0,k-deg2>::value+1>
                 (
                  coef1+max<0,k-deg2>::value
                  )
                 .array()
                 *
                 vecmap<min<k,deg1>::value-max<0,k-deg2>::value+1>
                 (
                  coef2+k-min<k,deg1>::value
                  ).reverse().array()
                 ).sum();
                _PolyConv<deg1,deg2,k-1>::compute( coefout, coef1, coef2 );
            }
        };
        
        template <int deg1,int deg2>
        struct
        _PolyConv<deg1,deg2,-1>
        {
            static void
            compute( double *coefout, const double *coef1, const double *coef2 )
            {
            }
        };
        
        template <int deg1,int deg2>
        struct
        PolyConv
        {
            static void
            compute( Eigen::Matrix<double,deg1+deg2+1,1> &coefout, const Eigen::Matrix<double,deg1+1,1> &coef1, const Eigen::Matrix<double,deg2+1,1> &coef2 )
            {
                _PolyConv<deg1,deg2,deg1+deg2>::compute(coefout.data(),coef1.data(),coef2.data());
            }
        };

        template <>
        struct
        PolyConv<Eigen::Dynamic,Eigen::Dynamic>
        {
            static void
            compute( Eigen::VectorXd &coefout, const Eigen::VectorXd &coef1, const Eigen::VectorXd &coef2 )
            {
                int deg1 = coef1.rows()-1;
                int deg2 = coef2.rows()-1;
                
                for ( int k = deg1+deg2; k >= 0; k-- )
                {
                    coefout[k] =
                    (
                     coef1.block(
                                 std::max(0,k-deg2),0,
                                 std::min(k,deg1)-std::max(0,k-deg2)+1,1
                                 ).array()
                     *
                     coef2.block(
                                 k-std::min(k,deg1),0,
                                 std::min(k,deg1)-std::max(0,k-deg2)+1,1
                                 ).reverse().array()
                     ).sum();
                }
            }
        };

        template<int deg,int k>
        struct
        _PolyVal
        {
            static double
            compute( const Eigen::Matrix<double,deg+1,1> &coef, const double &x )
            {
                return coef[k] + x*_PolyVal<deg,k-1>::compute(coef,x);
            }
        };
        
        template<int deg>
        struct
        _PolyVal<deg,0>
        {
            static double
            compute( const Eigen::Matrix<double,deg+1,1> &coef, const double &x )
            {
                return coef[0];
            }
        };
        
        template<int deg>
        struct
        PolyVal
        {
            static double
            compute( const Eigen::Matrix<double,deg+1,1> &coef, const double &x )
            {
                return Internal::_PolyVal<deg,deg>::compute( coef, x );
            }
        };

        template<>
        struct
        PolyVal<Eigen::Dynamic>
        {
            static double
            compute( const Eigen::VectorXd &coef, const double &x )
            {
                int deg = coef.rows()-1;
                double val = 0;
                for ( int i = 0; i <= deg; i++ ) val = coef(i) + x*val;
                return val;
            }
        };

        template<int deg>
        struct BisectionMethod
        {
            static double compute(const Eigen::Matrix<double,deg+1,1> &coef, const double lower, const double upper )
            {
                double a = lower;
                double fa = PolyVal<deg>::compute(coef,a);
                char fapos = (fa>0);
                double b = upper;
                double c = 0.5*(a+b);
                
                for ( int i = 0; i < 24; i++ )
                {
                    const double fc = PolyVal<deg>::compute(coef,c);
                    const char fcpos = (fc>0);
                    if ( fcpos == fapos )
                    {
                        a = c;
                        fa = fc;
                    } else {
                        b = c;
                    }
                    c = 0.5*(a+b);
                }
                
                return c;
            }
        };

        template<>
        struct BisectionMethod<Eigen::Dynamic>
        {
            static double compute(const Eigen::VectorXd &coef, const double lower, const double upper )
            {
                double a = lower;
                double fa = PolyVal<Eigen::Dynamic>::compute(coef,a);
                char fapos = (fa>0);
                double b = upper;
                double c = 0.5*(a+b);
                
                for ( int i = 0; i < 24; i++ )
                {
                    const double fc = PolyVal<Eigen::Dynamic>::compute(coef,c);
                    const char fcpos = (fc>0);
                    if ( fcpos == fapos )
                    {
                        a = c;
                        fa = fc;
                    } else {
                        b = c;
                    }
                    c = 0.5*(a+b);
                }
                
                return c;
            }
        };
        
        template<int deg>
        struct SturmChain
        {
            Eigen::Matrix<double,deg+1,1> coef;
            SturmChain<deg-1> next;
        };
        
        template<>
        struct SturmChain<0>
        {
            Eigen::Matrix<double,1,1> coef;
        };
        
        template<int degin,int degout,int k>
        struct _GetSturmChain
        {
            static SturmChain<degout>& compute( SturmChain<degin> &sc )
            {
                return _GetSturmChain<degin-1,degout,k-1>::compute(sc.next);
            }
        };
        
        template<int degin,int degout>
        struct _GetSturmChain<degin,degout,0>
        {
            static SturmChain<degout>& compute( SturmChain<degin> &sc )
            {
                return sc;
            }
        };
        
        template<int degin,int degout>
        struct GetSturmChain
        {
            static SturmChain<degout>& compute( SturmChain<degin> &sc )
            {
                return _GetSturmChain<degin,degout,degin-degout>::compute(sc);
            }
        };
        
        template <int deg,int k>
        struct _ModPoly
        {
            static const void compute( SturmChain<deg> &sc )
            {
                const Eigen::Matrix<double,k+3,1> &u( GetSturmChain<deg,k+2>::compute(sc).coef );
                const Eigen::Matrix<double,k+2,1> &v( GetSturmChain<deg,k+1>::compute(sc).coef );
                Eigen::Matrix<double,k+3,1> r = u;
                const double s = v(0)/fabs(v(0));
                r(k+2) *= s;
                const int blocklength = k+1;
                r.block(1,0,blocklength,1) = s * r.block(1,0,blocklength,1) - r(0) * v.block(1,0,blocklength,1);
                r.block(2,0,blocklength,1) = s * r.block(2,0,blocklength,1) - r(1) * v.block(1,0,blocklength,1);
                Eigen::Matrix<double,k+1,1> &rout( GetSturmChain<deg,k>::compute(sc).coef );
                rout = r.block(2,0,k+1,1);
                double f = -fabs(rout(0));
                rout /= f;
                _ModPoly<deg,k-1>::compute(sc);
            }
        };
        
        template <int deg>
        struct _ModPoly<deg,0>
        {
            static const void compute( SturmChain<deg> &sc )
            {
                const Eigen::Matrix<double,3,1> &u( GetSturmChain<deg,2>::compute(sc).coef );
                const Eigen::Matrix<double,2,1> &v( GetSturmChain<deg,1>::compute(sc).coef );
                Eigen::Matrix<double,3,1> r = u;
                const double s = v(0)/fabs(v(0));
                r(2) *= s;
                r(1) = s * r(1) - r(0) * v(1);
                r(2) = s * r(2) - r(1) * v(1);
                GetSturmChain<deg,0>::compute(sc).coef[0] = -r(2);
            }
        };
        
        template <int deg>
        struct ModPoly
        {
            static const void compute( SturmChain<deg> &sc )
            {
                _ModPoly<deg,deg-2>::compute(sc);
            }
        };
        
        template <int deg>
        struct _SignChanges
        {
            static const void compute( SturmChain<deg> &sc, const double x, const double lf, int &count )
            {
                double f = PolyVal<deg>::compute( sc.coef, x );
                count += (f*lf<0);
                _SignChanges<deg-1>::compute( sc.next, x, f, count );
            }
        };
        
        template <>
        struct _SignChanges<0>
        {
            static const void compute( SturmChain<0> &sc, const double x, const double lf, int &count )
            {
                double f = sc.coef[0];
                count += (f*lf<0);
            }
        };
        
        template <int deg>
        struct SignChanges
        {
            static const void compute( SturmChain<deg> &sc, const double x, int &count )
            {
                double f = PolyVal<deg>::compute( sc.coef, x );
                count = 0;
                _SignChanges<deg-1>::compute( sc.next, x, f, count );
            }
        };
        
        template <>
        struct SignChanges<0>
        {
            static const void compute( SturmChain<0> &sc, const double x, int &count )
            {
                count = 0;
            }
        };
        
        struct SturmInterval
        {
            double lb, ub;
            int sclb, scub;
            SturmInterval(const double _lb, const double _ub, const int _sclb, const int _scub)
            : lb(_lb), ub(_ub), sclb(_sclb), scub(_scub) { }
        };
        
        struct RootInterval
        {
            double lb, ub;
            RootInterval(const double _lb, const double _ub)
            : lb(_lb), ub(_ub) { }
        };
        
        template<int deg>
        class SturmRootFinder
        {
            void bisection(const double lb, const double ub, std::vector<double> &roots)
            {
                int sclb = 0;
                int scub = 0;
                SignChanges<deg>::compute(sc,lb,sclb);
                SignChanges<deg>::compute(sc,ub,scub);
                std::queue<SturmInterval> intervals;
                intervals.push( SturmInterval(lb,ub,sclb,scub) );
                
                std::vector<RootInterval> root_intervals;
                
                while ( !intervals.empty() )
                {
                    SturmInterval interval = intervals.front();
                    intervals.pop();
                    
                    // get num roots in interval
                    int n = interval.sclb-interval.scub;
                    
                    if ( n == 0 )
                    {
                        
                    }
                    else if ( n == 1 )
                    {
                        root_intervals.push_back( RootInterval(interval.lb,interval.ub) );
                    }
                    else if ( fabs(interval.ub-interval.lb)>1e-15 )
                    {
                        double m = (interval.lb+interval.ub)*0.5;
                        
                        // get sign changes at m
                        int scm = 0;
                        SignChanges<deg>::compute(sc,m,scm);
                        
                        // add [lb,m] interval
                        intervals.push( SturmInterval(interval.lb,m,interval.sclb,scm) );
                        
                        // add [m,ub] interval
                        intervals.push( SturmInterval(m,interval.ub,scm,interval.scub) );
                    }
                }
                
                roots.resize(root_intervals.size());
                for ( int i = 0; i < root_intervals.size(); i++ )
                {
                    roots[i] = BisectionMethod<deg>::compute(sc.coef,root_intervals[i].lb,root_intervals[i].ub);
                }
            }
        public:
            SturmChain<deg> sc;
            
            SturmRootFinder( const Eigen::Matrix<double,deg+1,1> &coef )
            {
                sc.coef = coef;
                sc.next.coef = coef.head(deg);
                double f = fabs(sc.coef[0] * deg);
                for ( int i = 0; i < deg; i++ ) sc.next.coef[i] *= (deg-i)/f;
                ModPoly<deg>::compute(sc);
            }
            
            void realRoots(const double lb, const double ub, std::vector<double> &roots)
            {
                bisection(lb,ub,roots);
            }
        };

        template<>
        class SturmRootFinder<Eigen::Dynamic>
        {
            std::vector<Eigen::VectorXd> sc;

            Eigen::VectorXd modpoly( const Eigen::VectorXd &u, const Eigen::VectorXd &v )
            {
                int udeg = u.rows()-1;
                Eigen::VectorXd r = u;
                double s = v(0)/fabs(v(0));
                r(udeg) *= s;
                int blocklength = udeg-1;
                r.block(1,0,blocklength,1) = s * r.block(1,0,blocklength,1) - r(0) * v.block(1,0,blocklength,1);
                r.block(2,0,blocklength,1) = s * r.block(2,0,blocklength,1) - r(1) * v.block(1,0,blocklength,1);
                return r.block(2,0,blocklength,1);
            }
            
            int signchanges(double x)
            {
                int deg = sc[0].rows()-1;
                Eigen::VectorXd f(deg+1);
                for ( int i = 0; i < deg+1; i++ ) f[i] = PolyVal<Eigen::Dynamic>::compute(sc[i],x);
                int n = 0;
                // NB: not skipping zeros here
                for ( int i = 1; i < deg+1; i++ ) n += (f[i]*f[i-1]<0);
                return n;
            }
            
            void bisection(const double lb, const double ub, std::vector<double> &roots)
            {
                int sclb = signchanges(lb);
                int scub = signchanges(ub);
                std::queue<SturmInterval> intervals;
                intervals.push( SturmInterval(lb,ub,sclb,scub) );
                
                std::vector<RootInterval> root_intervals;
                
                while ( !intervals.empty() )
                {
                    SturmInterval interval = intervals.front();
                    intervals.pop();
                    
                    // get num roots in interval
                    int n = interval.sclb-interval.scub;
                    
                    if ( n == 0 )
                    {
                        
                    }
                    else if ( n == 1 )
                    {
                        root_intervals.push_back( RootInterval(interval.lb,interval.ub) );
                    }
                    else if ( fabs(interval.ub-interval.lb)>1e-15 )
                    {
                        double m = (interval.lb+interval.ub)*0.5;
                        
                        // get sign changes at m
                        int scm = signchanges(m);
                        
                        // add [lb,m] interval
                        intervals.push( SturmInterval(interval.lb,m,interval.sclb,scm) );
                        
                        // add [m,ub] interval
                        intervals.push( SturmInterval(m,interval.ub,scm,interval.scub) );
                    }
                }
                
                roots.resize(root_intervals.size());
                for ( int i = 0; i < root_intervals.size(); i++ )
                {
                    roots[i] = BisectionMethod<Eigen::Dynamic>::compute(sc[0],root_intervals[i].lb,root_intervals[i].ub);
                }
            }
        public:
            SturmRootFinder( const Eigen::VectorXd &coef )
            {
                int deg = coef.rows()-1;
                sc.resize(deg+1);
                
                sc[0] = coef;
                
                double f = fabs(sc[0](0) * deg);
                sc[1] = sc[0].block(0,0,deg,1);
                for ( int i = 0; i < deg; i++ ) sc[1][i] *= (deg-i)/f;
                
                for ( int i = 2; i < deg+1; i++ )
                {
                    sc[i] = modpoly(sc[i-2],sc[i-1]);
                    if ( i != deg )
                    {
                        double f = -fabs(sc[i](0));
                        sc[i] /= f;
                    }
                    else
                    {
                        sc[i] = -sc[i];
                    }
                }
            }
            
            void realRoots(const double lb, const double ub, std::vector<double> &roots)
            {
                bisection(lb,ub,roots);
            }
        };

        template<int deg>
        struct RootFinder
        {
            static void compute( const Eigen::Matrix<double,deg+1,1> &p, std::vector<double> &roots )
            {
                Eigen::PolynomialSolver<double,deg> ps;
                ps.compute(p.reverse());
                ps.realRoots(roots);
            }
        };

        template<>
        struct RootFinder<Eigen::Dynamic>
        {
            static void compute( const Eigen::Matrix<double,Eigen::Dynamic,1> &p, std::vector<double> &roots )
            {
                Eigen::PolynomialSolver<double,Eigen::Dynamic> ps;
                ps.compute(p.reverse());
                ps.realRoots(roots);
            }
        };
        
        template<>
        struct RootFinder<2>
        {
            static void compute( const Eigen::Matrix<double,3,1> &p, std::vector<double> &roots )
            {
                const double discrim = p[1]*p[1] - 4.*p[0]*p[2];
                if ( discrim < 0 )
                {
                    roots.clear();
                }
                else if ( discrim == 0 )
                {
                    roots.resize(1);
                    roots[0] = -p[1] / ( 2. * p[0] );
                }
                else
                {
                    const double twoa = 2. * p[0];
                    const double sqrt_discrim = sqrt(discrim);
                    roots.resize(2);
                    roots[0] = (-p[1] + sqrt_discrim) / twoa;
                    roots[1] = (-p[1] - sqrt_discrim) / twoa;
                }
            }
        };

    } // end namespace Internal
    
} // end namespace Polynomial

#endif
