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
#ifndef itkIsSame_h
#define itkIsSame_h

namespace itk
{
  /** \cond HIDE_META_PROGRAMMING */
  /** True Type.
   * borrowed from type_traits
   * \ingroup MPL
   */
  struct TrueType
  {
    typedef bool     ValueType;
    typedef TrueType Type;

    static const ValueType Value = true;
    operator ValueType() { return Value; }
  };

  /** False type.
   * borrowed from type_traits
   * \ingroup MPL
   */
  struct FalseType
  {
    typedef bool      ValueType;
    typedef FalseType Type;
    static const ValueType Value = false;
    operator ValueType() { return Value; }
  };

  /** Tells whether two types are identical.
   * Default case: returns that types differ.
   * \ingroup MPL
   */
  template<typename, typename>
  struct IsSame
    : public FalseType
  {
  };

  /** \cond SPECIALIZATION_IMPLEMENTATION */
  /** Tells whether two types are identical.
   * Overload: returns that the two types are identical
   * \ingroup MPL
   */
  template<typename T>
  struct IsSame<T, T>
    : public TrueType
  {
  };
  /**\endcond*/

  /** \endcond */

} // end namespace itk

#endif //itkIsSame_h
