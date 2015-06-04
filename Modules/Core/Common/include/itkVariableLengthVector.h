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
#ifndef itkVariableLengthVector_h
#define itkVariableLengthVector_h

#define ITK_USE_EXPRESSION_TEMPLATE

#include <cassert>
#include "itkNumericTraits.h"

#if defined(ITK_USE_EXPRESSION_TEMPLATE)
#  include "itkPromoteType.h"
#  include "itkMPL.h"
#endif


namespace itk
{
  struct AlwaysReallocate {
      bool operator()(unsigned int newSize, unsigned int oldSize) const
      {
          (void)newSize;
          (void)oldSize;
          return true;
      }
  };
  struct NeverReallocate {
      bool operator()(unsigned int newSize, unsigned int oldSize) const
      {
          (void)newSize;
          (void)oldSize;
          assert(newSize == oldSize && "SetSize is expected to never change the VariableLengthVector size...");
          // A typical scenario would be:
          // VariableLengthVector<...> v;
          // v.SetSize(someFixedSize)
          // for (pixel : image) {
          //     assert(expression.size() == someFixedSize);
          //     v.FastAssign( expression );
          // }
          return true;
      }
  };
  struct ShrinkToFit {
      bool operator()(unsigned int newSize, unsigned int oldSize) const
      { return newSize != oldSize; }
  };
  struct DontShrinkToFit {
      bool operator()(unsigned int newSize, unsigned int oldSize) const
      { return newSize > oldSize; }
  };

  struct KeepOldValues {
    template <typename TValue>
      void operator()(unsigned int newSize, unsigned int oldSize, TValue * oldBuffer, TValue * newBuffer) const
      {
          const std::size_t nb = std::min(newSize, oldSize);
          std::copy(oldBuffer, oldBuffer+nb, newBuffer);
      }
  };
  struct DumpOldValues {
    template <typename TValue>
      void operator()(unsigned int newSize, unsigned int oldSize, TValue * oldBuffer, TValue * newBuffer) const
      {
          (void)oldSize;
          (void)newSize;
          (void)oldBuffer;
          (void)newBuffer;
      }
  };


#if defined(ITK_USE_EXPRESSION_TEMPLATE)
    template <typename TExpr1, typename TExpr2, typename  TBinaryOp>
        struct VLVExpr;
#endif

/** \class VariableLengthVector
 * \brief Represents an array whose length can be defined at run-time.
 *
 * This class is templated over the data type. This data-type is meant
 * to be a scalar, such as float, double etc...
 *
 * \note
 * ITK itself provides several classes that can serve as \c Arrays.
 * \li FixedArray - Compile time fixed length arrays that's intended to
 * represent an enumerated collection of \c n entities.
 *
 * \li Array - Run time resizeable array that is intended to hold a
 * collection of \c n entities
 *
 * \li Vector - Compile time fixed length array that is intended to hold
 * a collection of \c n data types. A vector usually has a mathematical meaning.
 * It should only be used when mathematical operations such as addition,
 * multiplication by a scalar, product etc make sense.
 *
 * \li VariableLengthVector - Run time array that is intended to hold a collection
 * of scalar data types. Again, it should be used only when mathematical
 * operations on it are relevant. If not, use an Array.
 *
 * \li Point - Represents the spatial coordinates of a spatial location. Operators
 * on Point reflect geometrical concepts.
 *
 * \par For the reasons listed above, you cannot instantiate
 * \code VariableLengthVector< bool > \endcode.
 *
 * \par
 * Design Considerations: We do not derive from \c vnl_vector to avoid being
 * limited by the explicit template instantiations of vnl_vector and other
 * hacks that vnl folks have been forced to use.
 *
 * \note
 * This work is part of the National Alliance for Medical Image Computing
 * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
 * for Medical Research, Grant U54 EB005149.
 *
 * \sa CovariantVector
 * \sa SymmetricSecondRankTensor
 * \sa RGBPixel
 * \sa DiffusionTensor3D
 * \ingroup DataRepresentation
 * \ingroup ITKCommon
 *
 * \wiki
 * \wikiexample{SimpleOperations/VariableLengthVector,Variable length vector}
 * \endwiki
 */
template< typename TValue >
class VariableLengthVector
{
public:

  /** The element type stored at each location in the Array. */
  typedef TValue                                        ValueType;
  typedef TValue                                        ComponentType;
  typedef typename NumericTraits< ValueType >::RealType RealValueType;
  typedef VariableLengthVector                          Self;

  /** Typedef used to indicate the number of elements in the vector */
  typedef unsigned int ElementIdentifier;

  /** Default constructor. It is created with an empty array
   *  it has to be allocated later by assignment              */
  VariableLengthVector();

  /** Constructor with size. Size can only be changed by assignment */
  explicit VariableLengthVector(unsigned int dimension);

  /** Constructor that initializes array with contents from a user supplied
   * buffer. The pointer to the buffer and the length is specified. By default,
   * the array does not manage the memory of the buffer. It merely points to
   * that location and it is the user's responsibility to delete it.
   * If "LetArrayManageMemory" is true, then this class will free the
   * memory when this object is destroyed. */
  VariableLengthVector(ValueType *data, unsigned int sz,
                       bool LetArrayManageMemory = false);

  /** Constructor that initializes array with contents from a user supplied
   * buffer. The pointer to the buffer and the length is specified. By default,
   * the array does not manage the memory of the buffer. It merely points to
   * that location and it is the user's responsibility to delete it.
   * If "LetArrayManageMemory" is true, then this class will free the
   * memory when this object is destroyed. */
  VariableLengthVector(const ValueType *data, unsigned int sz,
                       bool LetArrayManageMemory = false);

  /** Copy constructor. The reason why the copy constructor and the assignment
   * operator are templated is that it will allow implicit casts to be
   * performed. For instance
   * \code
   * VariableLengthVector< int > vI;
   * VariableLengthVector< float > vF( vI );
   * or for instance vF = static_cast< VariableLengthVector< float > >( vI );
   * \endcode
   */
  template< typename T >
  VariableLengthVector(const VariableLengthVector< T > & v)
  {
    m_NumElements = v.Size();
    m_Data = this->AllocateElements(m_NumElements);
    m_LetArrayManageMemory = true;
    for ( ElementIdentifier i = 0; i < m_NumElements; ++i )
      {
      this->m_Data[i] = static_cast< ValueType >( v[i] );
      }
  }

  /** Copy constructer.. Override the default non-templated copy constructor
   * that the compiler provides */
  VariableLengthVector(const VariableLengthVector< TValue > & v);

  /** Swaps two \c VariableLengthVector 's.
   * \pre Expects none of the \c VariableLengthVector to act as a proxy,
   * checked with assertions.
   * \param[in,out] v  other \c VariableLengthVector to be swapped with.
   * \throw None
   */
  void Swap(Self & v) {
      assert(m_LetArrayManageMemory);
      assert(v.m_LetArrayManageMemory);
      using std::swap;
      swap(v.m_Data       , m_Data);
      swap(v.m_NumElements, m_NumElements);
  }

#if defined(ITK_USE_EXPRESSION_TEMPLATE)
  template <typename TExpr1, typename TExpr2, typename  TBinaryOp>
      VariableLengthVector(VLVExpr<TExpr1, TExpr2, TBinaryOp> const& rhs_);
  template <typename TExpr1, typename TExpr2, typename  TBinaryOp>
  Self & operator=(const VLVExpr<TExpr1, TExpr2, TBinaryOp> & v);
#endif

  /** Set the all the elements of the array to the specified value */
  void Fill(TValue const & v);

  /** Assignment operator.
   * Oddly enough when the VLV owns the memory, old elements are not copied
   * ....
   */
  template< typename T >
  Self & operator=(const VariableLengthVector< T > & v)
  {
#if 0
      if ( m_Data == static_cast< void * >( const_cast< T * >
                  ( ( const_cast< VariableLengthVector<T> & >( v ) ).GetDataPointer() ) ) )
      {
          return *this;
      }
#endif
    ElementIdentifier const N = v.Size();
    this->SetSize( N, DontShrinkToFit(), DumpOldValues() );
    for ( ElementIdentifier i = 0; i < N; ++i )
      {
      this->m_Data[i] = static_cast< ValueType >( v[i] );
      }
    return *this;
  }

  /** Assignment operators  */
  Self & operator=(const Self & v);

  /** Fast Assignment */
  Self & FastAssign(const Self & v);

  Self & operator=(TValue const & v);

  /** Return the number of elements in the Array  */
  inline unsigned int Size(void) const { return m_NumElements; }
  inline unsigned int GetNumberOfElements(void) const { return m_NumElements; }

  /** Return reference to the element at specified index. No range checking. */
  TValue       & operator[](unsigned int i)       { return this->m_Data[i]; }
  /** Return reference to the element at specified index. No range checking. */
  TValue const & operator[](unsigned int i) const { return this->m_Data[i]; }

  /** Get one element */
  inline const TValue & GetElement(unsigned int i) const { return m_Data[i]; }

  /** Set one element */
  void SetElement(unsigned int i, const TValue & value) { m_Data[i] = value; }

  template <typename TReallocatePolicy, typename TKeepValuesPolicy>
  void SetSize(unsigned int sz,
          TReallocatePolicy reallocatePolicy,
          TKeepValuesPolicy keepValues);

  /** Set the size to that given.
   *
   * If \c destroyExistingData is \c false:
   * If the array already contains data, the existing data is copied over and
   * new space is allocated, if necessary. If the length to reserve is less
   * than the current number of elements, then an appropriate number of elements
   * are discarded.
   *    If \c true, the size is set destructively to the length given. If the
   * length is different from the current length, existing data will be lost.
   * The default is \c true. */
  void SetSize(unsigned int sz, bool destroyExistingData = true)
  {
      // Stays compatiable with previous code version
      // And works around the fact template function can't have default
      // arguments on template types.
      if (destroyExistingData)
          SetSize(sz, AlwaysReallocate(), KeepOldValues());
      else
          SetSize(sz, ShrinkToFit(), KeepOldValues());
  }

  /** Destroy data that is allocated internally, if LetArrayManageMemory is
   * true. */
  void DestroyExistingData();

  inline unsigned int GetSize(void) const { return m_NumElements; }

  /** Set the pointer from which the data is imported.
   * If "LetArrayManageMemory" is false, then the application retains
   * the responsibility of freeing the memory for this data.  If
   * "LetArrayManageMemory" is true, then this class will free the
   * memory when this object is destroyed. */
  void SetData(TValue *data, bool LetArrayManageMemory = false);

  /** Similar to the previous method. In the above method, the size must be
   * separately set prior to using user-supplied data. This introduces an
   * unnecessary allocation step to be performed. This method avoids it
   * and should be used to import data wherever possible to avoid this.
   * Set the pointer from which the data is imported.
   * If "LetArrayManageMemory" is false, then the application retains
   * the responsibility of freeing the memory for this data.  If
   * "LetArrayManageMemory" is true, then this class will free the
   * memory when this object is destroyed. */
  void SetData(TValue *data, unsigned int sz, bool LetArrayManageMemory = false);

  /** This destructor is not virtual for performance reasons. However, this
   * means that subclasses cannot allocate memory. */
  ~VariableLengthVector();

  /** Reserves memory of a certain length.
   *
   * If the array already contains data, the existing data is copied over and
   * new space is allocated, if necessary. If the length to reserve is less
   * than the current number of elements, then an appropriate number of elements
   * are discarded. */
  void Reserve(ElementIdentifier);

  /** Allocate memory of certain size and return it.  */
  TValue * AllocateElements(ElementIdentifier size) const;

  const TValue * GetDataPointer() const { return m_Data; }

#if !defined(ITK_USE_EXPRESSION_TEMPLATE)
  /** Element-wise vector addition. The vectors do not have to have
   * the same element type. The input vector elements are cast to the
   * output vector element type before the addition is performed.
   *
   * \note For efficiency, the length of the vectors is not checked;
   * they are assumed to have the same length. */
  template< typename T >
      inline Self operator+(const VariableLengthVector< T > & v) const
      {
          assert( m_NumElements == v.GetSize() );
          const ElementIdentifier length = v.Size();
          Self                    result(length);

          for ( ElementIdentifier i = 0; i < length; i++ )
          {
              result[i] = ( *this )[i] + static_cast< ValueType >( v[i] );
          }
          return result;
      }

  /** Element-wise subtraction of vectors. The vectors do not have to
   * have the same element type. The input vector elements are cast to
   * the output vector element type before the subtraction is
   * performed.
   *
   * \note For efficiency, the length of the vectors is not checked;
   * they are assumed to have the same length. */
  template< typename T >
  inline Self operator-(const VariableLengthVector< T > & v) const
  {
    assert( m_NumElements == v.GetSize() );
    const ElementIdentifier length = v.Size();
    Self                    result(length);

    for ( ElementIdentifier i = 0; i < length; i++ )
      {
      result[i] = ( *this )[i] - static_cast< ValueType >( v[i] );
      }
    return result;
  }

  /** Multiply vector elements by a scalar 's'. The vector does not
   * have to have the same element type as the scalar type. The scalar
   * is cast to the output vector element type before the
   * multiplication is performed. */
  template< typename T >
  inline Self operator*(T s) const
  {
    Self result(m_NumElements);

    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      result[i] = m_Data[i] * static_cast< ValueType >( s );
      }
    return result;
  }

  /** Divide vector elements by a scalar 's'. The vector does not
   * have to have the same element type as the scalar type. Both the
   * scalar and vector elements are cast to the RealValueType prior to
   * division, and the result is cast to the ValueType. */
  template< typename T >
  inline Self operator/(T s) const
  {
    Self result(m_NumElements);

    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      result[i] = static_cast< ValueType >(
        static_cast< RealValueType >( m_Data[i] )
        / static_cast< RealValueType >( s ) );
      }
    return result;
  }

  /** Add scalar 's' to each element of the vector.*/
  inline Self operator+(TValue s) const
  {
      Self result(m_NumElements);

      for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
          result[i] = m_Data[i] + s;
      }
      return result;
  }

  /** Subtract scalar 's' from each element of the vector.*/
  inline Self operator-(TValue s) const
  {
    Self result(m_NumElements);

    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      result[i] = m_Data[i] - s;
      }
    return result;
  }
#endif

  /** Prefix operator that subtracts 1 from each element of the
   * vector. */
  inline Self & operator--()
  {
    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      this->m_Data[i] -= static_cast< ValueType >( 1.0 );
      }
    return *this;
  }

  /** Prefix operator that adds 1 to each element of the vector. */
  inline Self & operator++() // prefix operator ++v;
  {
    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      this->m_Data[i] += static_cast< ValueType >( 1.0 );
      }
    return *this;
  }

  /** Postfix operator that subtracts 1 from each element of the
   * vector. */
  inline Self operator--(int) // postfix operator v--;
  {
    Self tmp(*this);

    --tmp;
    return tmp;
  }

  /** Postfix operator that adds 1 to each element of the vector. */
  inline Self operator++(int) // postfix operator v++;
  {
    Self tmp(*this);

    ++tmp;
    return tmp;
  }

  /** Element-wise subtraction of vector 'v' from the current
   * vector. The vectors do not have to have the same element
   * type. The input vector elements are cast to the current vector
   * element type before the subtraction is performed.
   *
   * \note For efficiency, the length of the vectors is not checked;
   * they are assumed to have the same length. */
  template< typename T >
  inline Self & operator-=
    (const VariableLengthVector< T > & v)
  {
    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      m_Data[i] -= static_cast< ValueType >( v[i] );
      }
    return *this;
  }

  /** Subtract scalar 's' from each element of the current vector. */
  inline Self & operator-=(TValue s)
  {
    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      m_Data[i] -= s;
      }
    return *this;
  }

  /** Element-wise addition of vector 'v' to the current vector. The
   * vectors do not have to have the same element type. The input
   * vector elements are cast to the current vector element type
   * before the addition is performed.
   *
   * \note For efficiency, the length of the vectors is not checked;
   * they are assumed to have the same length. */
  template< typename T >
  inline Self & operator+=
    (const VariableLengthVector< T > & v)
  {
    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      m_Data[i] += static_cast< ValueType >( v[i] );
      }
    return *this;
  }

  /** Add scalar 's' to each element of the vector. */
  inline Self & operator+=(TValue s)
  {
    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      m_Data[i] += s;
      }
    return *this;
  }

#if defined(ITK_USE_EXPRESSION_TEMPLATE)
  template <typename TExpr1, typename TExpr2, typename TBinaryOp>
    Self& operator+=(VLVExpr<TExpr1,TExpr2,TBinaryOp> const& v)
      {
      for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
        {
        m_Data[i] += static_cast< ValueType >( v[i] );
        }
      return *this;
      }
#endif

  /** Multiply each element of the vector by a scalar 's'. The scalar
   * value is cast to the current vector element type prior to
   * multiplication. */
  template< typename T >
  inline Self & operator*=(T s)
  {
    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      m_Data[i] *= ( static_cast< ValueType >( s ) );
      }
    return *this;
  }

  /** Divide vector elements by a scalar 's'. The vector does not
   * have to have the same element type as the scalar type. Both the
   * scalar and vector elements are cast to the RealValueType prior to
   * division, and the result is cast to the ValueType. */
  template< typename T >
  inline Self & operator/=(T s)
  {
    for ( ElementIdentifier i = 0; i < m_NumElements; i++ )
      {
      m_Data[i] = static_cast< ValueType >(
        static_cast< RealValueType >( m_Data[i] )
        / static_cast< RealValueType >( s ) );
      }
    return *this;
  }

  /** Negates each vector element. */
  Self & operator-();  // negation operator

  bool operator==(const Self & v) const;

  bool operator!=(const Self & v) const;

  /** Returns vector's Euclidean Norm  */
  RealValueType GetNorm() const;

  /** Returns vector's squared Euclidean Norm  */
  RealValueType GetSquaredNorm() const;

  /** letArrayManageMemory getter. */
  bool isAProxy() const { return ! m_LetArrayManageMemory;}
private:

  bool              m_LetArrayManageMemory; // if true, the array is responsible
                                            // for memory of data
  TValue *          m_Data;                 // Array to hold data
  ElementIdentifier m_NumElements;
};

#if !defined(ITK_USE_EXPRESSION_TEMPLATE)
/** Premultiply Operator for product of a VariableLengthVector and a scalar.
 *  VariableLengthVector< TValue >  =  T * VariableLengthVector< TValue >
 */
template< typename TValue, typename T >
inline
VariableLengthVector< TValue >
operator*(const T & scalar, const VariableLengthVector< TValue > & v)
{
  return v.operator*(scalar);
}
#endif

template< typename TValue >
std::ostream & operator<<(std::ostream & os, const VariableLengthVector< TValue > & arr)
{
  const unsigned int length = arr.Size();
  const signed int   last   = (unsigned int)length - 1;

  os << "[";
  for ( signed int i = 0; i < last; ++i )
    {
    os << arr[i] << ", ";
    }
  if ( length >= 1 )
    {
    os << arr[last];
    }
  os << "]";
  return os;
}

#if defined(ITK_USE_EXPRESSION_TEMPLATE)
template <typename TExpr> struct get_type
{
    typedef TExpr Type;
    static type load(type const& v, unsigned int idx)
    { (void)idx; return v; }
};
template <typename TExpr1, typename TExpr2> struct get_size
{
    static unsigned int Size(TExpr1 const&, TExpr2 const& rhs_)
    { return rhs_.Size(); }
};
template <typename TExpr1, typename TExpr2, typename TBinaryOp, typename  TRHS>
struct get_size<VLVExpr<TExpr1, TExpr2, TBinaryOp>, TRHS>
{
    static unsigned int Size(VLVExpr<TExpr1, TExpr2, TBinaryOp> const& lhs_, TRHS const& )
    { return lhs_.Size(); }
};
template <typename T, typename  TRHS>
struct get_size<VariableLengthVector<T>, TRHS>
{
    static unsigned int Size(VariableLengthVector<T> const& lhs_, TRHS const& )
    { return lhs_.Size(); }
};

template <typename TExpr1, typename TExpr2, typename  TBinaryOp>
struct VLVExpr
{
    static const bool isVLV = true;

    VLVExpr(TExpr1 const& lhs_, TExpr2 const& rhs_) : m_lhs(lhs_), m_rhs(rhs_){}
    // TODO: detect constants to avoid computing sizes...
    unsigned int Size() const { return get_size<TExpr1, TExpr2>::Size(m_lhs, m_rhs); }

    typedef typename  PromoteType<
        typename get_type<TExpr1>::Type,
        typename get_type<TExpr2>::Type        >::type  ResType;
    typedef typename NumericTraits< ResType >::RealType RealValueType;

    ResType operator[](unsigned int idx) const {
        return TBinaryOp::apply(
                get_type<TExpr1>::load(m_lhs, idx),
                get_type<TExpr2>::load(m_rhs, idx));
    }

    /** Returns vector's Euclidean Norm  */
    RealValueType GetNorm() const;

    /** Returns vector's squared Euclidean Norm  */
    RealValueType GetSquaredNorm() const;

private:
    TExpr1 const& m_lhs;
    TExpr2 const& m_rhs;
};

template <typename T>
struct get_type<VariableLengthVector<T> >
{
    typedef T Type;
    static type load(VariableLengthVector<T>  const& v, unsigned int idx)
    { return v[idx]; }
};
template <typename TExpr1, typename TExpr2, typename TBinaryOp>
struct get_type<VLVExpr<TExpr1, TExpr2, TBinaryOp> >
{
    typedef typename VLVExpr<TExpr1, TExpr2, TBinaryOp>::ResType Type;
    static type load(VLVExpr<TExpr1, TExpr2, TBinaryOp> const& v, unsigned int idx)
    { return v[idx]; }
};

namespace op {
    struct plus  {
        template <typename T1, typename T2>
            static
            typename PromoteType<T1, T2>::type apply(T1 const& lhs_, T2 const& rhs_)
            { return lhs_ + rhs_; }
    };

    struct sub  {
        template <typename T1, typename T2>
            static
            typename PromoteType<T1, T2>::type apply(T1 const& lhs_, T2 const& rhs_)
            { return lhs_ - rhs_; }
    };

    struct mult  {
        template <typename T1, typename T2>
            static
            typename PromoteType<T1, T2>::type apply(T1 const& lhs_, T2 const& rhs_)
            { return lhs_ * rhs_; }
    };

    struct div  {
        template <typename T1, typename T2>
            static
            typename PromoteType<T1, T2>::type apply(T1 const& lhs_, T2 const& rhs_)
            { return lhs_ / rhs_; }
    };
} // op namespace

#if 1
template <typename T>
struct is_VLV : false_type {};

template <typename T>
struct is_VLV<itk::VariableLengthVector<T> > : true_type {};

template <typename TExpr1, typename TExpr2, typename TBinaryOp>
struct is_VLV<VLVExpr<TExpr1, TExpr2,TBinaryOp> > : true_type {};
#endif

template <typename TExpr1, typename TExpr2>
inline
typename enable_if<mpl::or_<is_VLV<TExpr1>, is_VLV<TExpr2> > , VLVExpr<TExpr1, TExpr2, op::plus> >::type
operator+(TExpr1 const& lhs_, TExpr2 const& rhs_)
{ return VLVExpr<TExpr1, TExpr2, op::plus>(lhs_, rhs_); }

template <typename TExpr1, typename TExpr2>
inline
typename enable_if<mpl::or_<is_VLV<TExpr1>, is_VLV<TExpr2> > , VLVExpr<TExpr1, TExpr2, op::sub> >::type
operator-(TExpr1 const& lhs_, TExpr2 const& rhs_)
{ return VLVExpr<TExpr1, TExpr2, op::sub>(lhs_, rhs_); }

template <typename TExpr1, typename TExpr2>
inline
typename enable_if<mpl::xor_<is_VLV<TExpr1>, is_VLV<TExpr2> > , VLVExpr<TExpr1, TExpr2, op::mult> >::type
operator*(TExpr1 const& lhs_, TExpr2 const& rhs_)
{ return VLVExpr<TExpr1, TExpr2, op::mult>(lhs_, rhs_); }

template <typename TExpr1, typename TExpr2>
inline
typename enable_if<mpl::and_<is_VLV<TExpr1>, mpl::not_<is_VLV<TExpr2> > > , VLVExpr<TExpr1, TExpr2, op::div> >::type
operator/(TExpr1 const& lhs_, TExpr2 const& rhs_)
{ return VLVExpr<TExpr1, TExpr2, op::div>(lhs_, rhs_); }

template <typename TExpr1, typename TExpr2, typename  TBinaryOp>
std::ostream & operator<<(std::ostream &os, VLVExpr<TExpr1, TExpr2, TBinaryOp> const& v)
{
    os << "[";
    if (v.Size() != 0) {
        os << v[0];
        for (unsigned int i = 1, N = v.Size(); i != N; ++i) {
            os << ", " << v[i];
        }
    }
    return os << "]";
}

template <typename TExpr>
inline
typename enable_if<is_VLV<TExpr>, typename TExpr::RealValueType>::type
GetNorm(TExpr const& v)
{ return static_cast<typename TExpr::RealValueType>(std::sqrt(static_cast<double>(GetSquaredNorm(v)))); }

template <typename TExpr>
inline
typename enable_if<is_VLV<TExpr>, typename TExpr::RealValueType>::type
GetSquaredNorm(TExpr const& v)
{
  typedef typename TExpr::RealValueType RealValueType;
  RealValueType sum = 0.0;
  for ( unsigned int i = 0, N=v.Size(); i < N; ++i )
    {
    const RealValueType value = v[i];
    sum += value * value;
    }
  return sum;
}
#endif

template <typename T>
void swap(VariableLengthVector<T> &l_, VariableLengthVector<T> &r_)
{ l_.Swap(r_); }

// MACCS 4.3 - optimization
// aim: avoid a static cast that generates a new VLV
template <typename TDest, typename TOrig>
inline
void CastInto(VariableLengthVector<TDest> & d, VariableLengthVector<TOrig> const& o)
{ d = o; }

/**
 * Specialized function to convert a vector pixel into another vector pixel type.
 * This function is aimed at writing efficient algorithms that move pixel
 * values around. It is particularly important to not degrade performances on
 * \c itk::VariableLentghVector based pixels.
 *
 * While this overload does not really move (the internals of) \c o into \d
 * (as C++11 rvalue-references would permit), it will nonetheless avoid
 * allocating memory for a new \c VariableLentghVector, thing that would have
 * been a consequence of a \c static_cast<>.
 *
 * Instead, this overload will use the conversion assignment operator from \c
 * VariableLentghVector to copy (and convert on-the-fly) the value from \c o.

 * @tparam TDest Internal pixel type of the destination pixel
 * @tparam TOrig Internal pixel type of the origin pixel
 * @param[out] d Destination pixel that shal receive the value of \c o
 * @param[in]  o Origin pixel
 * @pre \c d shall not be a \c VariableLentghVector that acts as a proxy.
 *
 * @throw None
 * \see <tt>itk::MoveInto(TDest, TOrig)</tt>
 * \see <tt>itk::MoveInto(VariableLengthVector<T>, VariableLengthVector<T>)</tt>
 * \ingroup DataRepresentation
 * \ingroup ITKCommon
 */
template <typename TDest, typename TOrig>
inline
void MoveInto(VariableLengthVector<TDest> & d, VariableLengthVector<TOrig> const& o)
{ d = o; }

// MACCS 4.3 - optimization
/**
 * Specialized function to move a vector pixel into another vector pixel.
 * This function is aimed at writing efficient algorithms that move pixel
 * values around. It is particularly important to not degrade performances on
 * \c itk::VariableLentghVector based pixels.
 *
 * Unlike the other overload that works on \c VariableLentghVector, this one
 * will move the content of \c o into \d. As we don't care about old \c o
 * content. Morevover in order to avoid the endless allocate-free cycles, the
 * old content of \c o is moved into \c d. In other word, \c o and \c d
 * contents are swapped.
 *
 * @tparam TDest Internal pixel type of the destination pixel
 * @tparam TOrig Internal pixel type of the origin pixel
 * @param[in,out] d Destination pixel that shal receive the value of \c o
 * @param[in,out] o Origin pixel, that will received old \c d value.
 * @pre Neither \c d nor \c o shall not be \c VariableLentghVector s that act as proxy.
 *
 * @throw None
 * \see <tt>itk::MoveInto(VariableLengthVector<TDest>, VariableLengthVector<TOrig>)</tt>
 * \see <tt>itk::MoveInto(TDest, TOrig)</tt>
 * \ingroup DataRepresentation
 * \ingroup ITKCommon
 */
template <typename T>
inline
void MoveInto(VariableLengthVector<T> & d, VariableLengthVector<T> & o)
{ swap(d, o); }

} // namespace itk

#include "itkNumericTraitsVariableLengthVectorPixel.h"

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVariableLengthVector.hxx"
#endif

#endif
