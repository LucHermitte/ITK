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

#include "itkIsSame.h"
#include "itkEnableIf.h"

namespace itk {
  /** \cond HIDE_META_PROGRAMMING */
  /**\defgroup MPL Meta Programming Library
   * \ingroup Common
   * This module contains definitions aimed at metaprogramming.
   * They are mainly codes borrowed or inspired from C++11 or indirectly boost,
   * with with ITK UpperCamelCase naming policy.
   */


  /** MPL \c if control-statement.
   * \tparam VP Boolean predicate
   * \tparam T1 Type returned if the predicate is true
   * \tparam T2 Type returned if the predicate is false
   * \ingroup MPL
   */
  template <bool VP, typename T1, typename T2> struct If_;
  /** \cond SPECIALIZATION_IMPLEMENTATION */
  template <typename T1, typename T2> struct If_<true , T1, T2>{ typedef T1 Type; };
  template <typename T1, typename T2> struct If_<false, T1, T2>{ typedef T2 Type; };
  /**\endcond*/

  namespace mpl {

    /** MPL \c OR operator on constants.
     * \tparam VF1 First boolean expression
     * \tparam VF2 Second boolean expression
     * \ingroup MPL
     */
    template < bool VF1, bool VF2> struct OrC : TrueType { };
    /** \cond SPECIALIZATION_IMPLEMENTATION */
    template <> struct OrC<false, false> : FalseType {};
    /**\endcond*/
    /** MPL \c OR operator on types.
     * \tparam TF1 First boolean type
     * \tparam TF2 Second boolean type
     *
     * This \em overload automatically fetches \c TF1 and \c TF2 values.
     * However, beware, it won't work with standard C++ traits or boost traits.
     * Indeed, this overload expects the \em value to follow UpperCamelCase ITK
     * naming policy instead of the standard snake_case policy.
     * \ingroup MPL
     */
    template < typename TF1, typename TF2> struct Or_ : OrC<TF1::Value, TF2::Value>
    { typedef typename OrC<TF1::Value, TF2::Value>::Type Type; };

    /** MPL \c AND operator on constants.
     * \tparam VF1 First boolean expression
     * \tparam VF2 Second boolean expression
     * \ingroup MPL
     */
    template < bool VF1, bool VF2> struct AndC : FalseType { };
    /** \cond SPECIALIZATION_IMPLEMENTATION */
    template <> struct AndC<true, true> : TrueType {};
    /**\endcond*/
    /** MPL \c AND operator on types.
     * \tparam TF1 First boolean type
     * \tparam TF2 Second boolean type
     *
     * This \em overload automatically fetches \c TF1 and \c TF2 values.
     * However, beware, it won't work with standard C++ traits or boost traits.
     * Indeed, this overload expects the \em value to follow UpperCamelCase ITK
     * naming policy instead of the standard snake_case policy.
     * \ingroup MPL
     */
    template < typename TF1, typename TF2> struct And_ : AndC<TF1::Value, TF2::Value>
    { typedef typename AndC<TF1::Value, TF2::Value>::Type Type; };

    /** MPL \c XOR operator on constants.
     * \tparam VF1 First boolean expression
     * \tparam VF2 Second boolean expression
     * \ingroup MPL
     */
    template < bool VF1, bool VF2> struct XorC : FalseType { };
    /** \cond SPECIALIZATION_IMPLEMENTATION */
    template <> struct XorC<true, false> : TrueType {};
    template <> struct XorC<false, true> : TrueType {};
    /**\endcond*/
    /** MPL \c XOR operator on types.
     * \tparam TF1 First boolean type
     * \tparam TF2 Second boolean type
     *
     * This \em overload automatically fetches \c TF1 and \c TF2 values.
     * However, beware, it won't work with standard C++ traits or boost traits.
     * Indeed, this overload expects the \em value to follow UpperCamelCase ITK
     * naming policy instead of the standard snake_case policy.
     * \ingroup MPL
     */
    template < typename TF1, typename TF2> struct Xor_ : XorC<TF1::Value, TF2::Value>
    { typedef typename XorC<TF1::Value, TF2::Value>::Type Type; };

    /** MPL \c NOT operator on constants.
     * \tparam VF boolean expression to negate
     * \ingroup MPL
     */
    template < bool VF > struct NotC : FalseType { };
    /** \cond SPECIALIZATION_IMPLEMENTATION */
    template <> struct NotC<false> : TrueType {};
    template <> struct NotC<true>  : FalseType {};
    /**\endcond*/
    /** MPL \c NOT operator on types.
     * \tparam TF Second boolean type
     *
     * This \em overload automatically fetches \c TF value.
     * However, beware, it won't work with standard C++ traits or boost traits.
     * Indeed, this overload expects the \em value to follow UpperCamelCase ITK
     * naming policy instead of the standard snake_case policy.
     * \ingroup MPL
     */
    template < typename TF> struct Not_ : NotC<TF::Value>
    { typedef typename NotC<TF::Value>::Type Type; };

  } // mpl namespace

    /** \endcond */
} // itk namespace
#endif // itkMPL_h
