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
        template <int N, typename TA, typename TB> struct size_to_type;
        template <int N> struct identity { typedef char type[N]; };

#define ASSOC(N, typed)\
        template <typename TA, typename TB> struct size_to_type<N,TA,TB> { typedef typed type; };

        ASSOC(1, A);
        ASSOC(2, B);

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
        template <typename T> static details::identity<1>::type& check(typename details::size_to_type<1,  TA, TB>::type, T);
        template <typename T> static details::identity<2>::type& check(typename details::size_to_type<2,  TA, TB>::type, T);

        // Common numeric types
        static details::identity<3 >::type& check(typename details::size_to_type<3,  TA, TB>::type, int);
        static details::identity<4 >::type& check(typename details::size_to_type<4,  TA, TB>::type, int);
        static details::identity<5 >::type& check(typename details::size_to_type<5,  TA, TB>::type, int);
        static details::identity<6 >::type& check(typename details::size_to_type<6,  TA, TB>::type, int);
        static details::identity<7 >::type& check(typename details::size_to_type<7,  TA, TB>::type, int);
        static details::identity<8 >::type& check(typename details::size_to_type<8,  TA, TB>::type, int);
        static details::identity<9 >::type& check(typename details::size_to_type<9,  TA, TB>::type, int);
        static details::identity<10>::type& check(typename details::size_to_type<10, TA, TB>::type, int);
        static details::identity<11>::type& check(typename details::size_to_type<11, TA, TB>::type, int);
        static details::identity<12>::type& check(typename details::size_to_type<12, TA, TB>::type, int);
        static details::identity<13>::type& check(typename details::size_to_type<13, TA, TB>::type, int);
    public:
        typedef typename details::size_to_type<sizeof check(a+b, 0), TA, TB>::type type;
    };


} // itk namespace

#if 0
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<int,int>         ::type, int>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<short,int>       ::type, int>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<double,int>      ::type, double>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<float,int>       ::type, float>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<long,int>        ::type, long>));
BOOST_MPL_ASSERT((boost::is_same<itk::PromoteType<long,long double>::type, long double>));
#endif

#endif // itkPromoteType_h
