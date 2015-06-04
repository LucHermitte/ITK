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
#ifndef itkMPL_h
#define itkMPL_h

// Excerpt of codes from
// - boost/utility/enable_if.hpp
// - boost/mpl/bool.hpp

namespace itk {
    template <bool TP, typename T1, typename T2> struct if_;
    template <typename T1, typename T2> struct if_<true , T1, T2>{ typedef T1 type; };
    template <typename T1, typename T2> struct if_<false, T1, T2>{ typedef T2 type; };

    // Copyright 2003 (c) The Trustees of Indiana University.

    // Use, modification, and distribution is subject to the Boost Software
    // License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    // http://www.boost.org/LICENSE_1_0.txt)

    //    Authors: Jaakko Jarvi (jajarvi at osl.iu.edu)
    //             Jeremiah Willcock (jewillco at osl.iu.edu)
    //             Andrew Lumsdaine (lums at osl.iu.edu)
    template <bool TB, class T = void>
        struct enable_if_c {
            typedef T type;
        };

    template <class T>
        struct enable_if_c<false, T> {};

    template <class TCond, class T = void>
        struct enable_if : public enable_if_c<TCond::value, T> {};

    //  (C) Copyright Steve Cleary, Beman Dawes, Howard Hinnant & John Maddock 2000.
    //  Use, modification and distribution are subject to the Boost Software License,
    //  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    //  http://www.boost.org/LICENSE_1_0.txt).
    //
    //  See http://www.boost.org/libs/type_traits for most recent version including documentation.

    template <class T, T Val>
        struct integral_constant
        {
            typedef integral_constant<T, Val> type;
            typedef T                         value_type;
            static const T value = Val;
        };

    typedef integral_constant<bool, true>  true_type;
    typedef integral_constant<bool, false> false_type;

    namespace mpl {

        template < bool TF1, bool TF2> struct or_c : true_type { };
        template <> struct or_c<false, false> : false_type {};
        template < typename TF1, typename TF2> struct or_ : or_c<TF1::value, TF2::value>
        { typedef typename or_c<TF1::value, TF2::value>::type type; };

        template < bool TF1, bool TF2> struct and_c : false_type { };
        template <> struct and_c<true, true> : true_type {};
        template < typename TF1, typename TF2> struct and_ : and_c<TF1::value, TF2::value>
        { typedef typename and_c<TF1::value, TF2::value>::type type; };

        template < bool TF1, bool TF2> struct xor_c : false_type { };
        template <> struct xor_c<true, false> : true_type {};
        template <> struct xor_c<false, true> : true_type {};
        template < typename TF1, typename TF2> struct xor_ : xor_c<TF1::value, TF2::value>
        { typedef typename xor_c<TF1::value, TF2::value>::type type; };

        template < bool TF > struct not_c : false_type { };
        template <> struct not_c<false> : true_type {};
        template <> struct not_c<true>  : false_type {};
        template < typename TF> struct not_ : not_c<TF::value>
        { typedef typename not_c<TF::value>::type type; };

    } // mpl namespace

} // itk namespace
#endif // itkMPL_h
