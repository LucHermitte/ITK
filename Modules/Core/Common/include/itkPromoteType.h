/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkPromoteType_h
#define itkPromoteType_h

// Simplification of boost::common_type
namespace itk {
    namespace details {
        template <int N, typename TA, typename TB> struct SizeToType;
        template <int N> struct Identity { typedef char Type[N]; };

#define ASSOC(N, Typed)\
        template <typename TA, typename TB> struct SizeToType<N,TA,TB> { typedef Typed Type; };

        ASSOC(1, TA);
        ASSOC(2, TB);

        ASSOC(3,  short);
        ASSOC(4,  unsigned short);
        ASSOC(5,  int);
        ASSOC(6,  unsigned int);
        ASSOC(7,  long);
        ASSOC(8,  unsigned long);
        ASSOC(9,  long long);
        ASSOC(10, unsigned long long);
        ASSOC(11, float);
        ASSOC(12, double);
        ASSOC(13, long double);
#undef ASSOC
    } // details namespace

    template <typename TA, typename TB> struct PromoteType
    {
        static TA a;
        static TB b;

        // Aimed at supporting overloads
        template <typename T> static details::Identity<1>::Type& check(typename details::SizeToType<1,  TA, TB>::Type, T);
        template <typename T> static details::Identity<2>::Type& check(typename details::SizeToType<2,  TA, TB>::Type, T);

        // Common numeric types
        static details::Identity<3 >::Type& check(typename details::SizeToType<3,  TA, TB>::Type, int);
        static details::Identity<4 >::Type& check(typename details::SizeToType<4,  TA, TB>::Type, int);
        static details::Identity<5 >::Type& check(typename details::SizeToType<5,  TA, TB>::Type, int);
        static details::Identity<6 >::Type& check(typename details::SizeToType<6,  TA, TB>::Type, int);
        static details::Identity<7 >::Type& check(typename details::SizeToType<7,  TA, TB>::Type, int);
        static details::Identity<8 >::Type& check(typename details::SizeToType<8,  TA, TB>::Type, int);
        static details::Identity<9 >::Type& check(typename details::SizeToType<9,  TA, TB>::Type, int);
        static details::Identity<10>::Type& check(typename details::SizeToType<10, TA, TB>::Type, int);
        static details::Identity<11>::Type& check(typename details::SizeToType<11, TA, TB>::Type, int);
        static details::Identity<12>::Type& check(typename details::SizeToType<12, TA, TB>::Type, int);
        static details::Identity<13>::Type& check(typename details::SizeToType<13, TA, TB>::Type, int);
    public:
        typedef typename details::SizeToType<sizeof check(a+b, 0), TA, TB>::Type Type;
    };


} // itk namespace

#if 0
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<int,int>         ::Type, int>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<short,int>       ::Type, int>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<double,int>      ::Type, double>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<float,int>       ::Type, float>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<long,int>        ::Type, long>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<long,long double>::Type, long double>));
#endif

#endif // itkPromoteType_h
