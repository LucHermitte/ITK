/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkRandomAccessNeighborhoodIterator.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) 2000 National Library of Medicine
  All rights reserved.

  See COPYRIGHT.txt for copyright details.

  =========================================================================*/
namespace itk {
  
template<class TPixel, unsigned int VDimension, class TAccessor,
    class TDerefAccessor>
RandomAccessNeighborhoodIterator<TPixel, VDimension, TAccessor, TDerefAccessor>
RandomAccessNeighborhoodIterator<TPixel, VDimension, TAccessor, TDerefAccessor>
::Begin() const
{
  //Copy the current iterator
  Self it( *this );

  // Set the position to the m_BeginOffset
  it.SetLocation( this->m_StartIndex );

  return it;
}

template<class TPixel, unsigned int VDimension, class TAccessor,
    class TDerefAccessor>
RandomAccessNeighborhoodIterator<TPixel, VDimension, TAccessor, TDerefAccessor>
RandomAccessNeighborhoodIterator<TPixel, VDimension, TAccessor, TDerefAccessor>
::End() const 
{
  IndexType endIndex;
  
  // Copy the current iterator
  Self it( *this );

  // Calculate the end index
  for (unsigned int i = 0; i< VDimension; ++i)
    {
      endIndex.m_Index[i] = m_Bound[i] -1;
    }
  it.SetLocation( endIndex );

  ++it;

  return it;
}


template<class TPixel, unsigned int VDimension, class TAccessor,
  class TDerefAccessor>
RandomAccessNeighborhoodIterator<TPixel, VDimension, TAccessor, TDerefAccessor> &
RandomAccessNeighborhoodIterator<TPixel, VDimension, TAccessor, TDerefAccessor>
::operator+=(const Index<VDimension> & idx)
{
  unsigned int i;
  Iterator it;
  const Iterator _end = this->end();
  unsigned long accumulator = 0;
  const unsigned long* stride = this->GetImagePointer()->GetOffsetTable();

  // Offset from the increment in the lowest dimension
  accumulator += idx[0];
  
  // Offsets from the stride lengths in each dimension.
  //
  // Because the image offset table is based on its buffer size and not its
  // requested region size, we don't have to worry about adding in the wrapping
  // offsets. 
  for (i = 1; i< VDimension; ++i)
    {
      accumulator += idx[i] * stride[i];
    }

  // Increment pointers.
  for (it = this->begin(); it < _end; ++it)
    {
      (*it) += accumulator;
    }
  if (m_OutputBuffer)
    {
      m_OutputBuffer += accumulator;
    }

  // Update loop counter values
  for (i=0; i<VDimension; ++i)
    {
      m_Loop[i]+= idx[i];
    }

  return *this;
}

template<class TPixel, unsigned int VDimension, class TAccessor,
  class TDerefAccessor>
RandomAccessNeighborhoodIterator<TPixel, VDimension, TAccessor,
  TDerefAccessor> &
RandomAccessNeighborhoodIterator<TPixel, VDimension, TAccessor, TDerefAccessor> 
::operator-=(const Index<VDimension> & idx)
{
  unsigned int i;
  Iterator it;
  const Iterator _end = this->end();
  unsigned long accumulator = 0;
  const unsigned long* stride = this->GetImagePointer()->GetOffsetTable();

  // Offset from the increment in the lowest dimension
  accumulator += idx[0];
  
  // Offsets from the stride lengths in each dimension.
  //
  // Because the image offset table is based on its buffer size and not its
  // requested region size, we don't have to worry about adding in the wrapping
  // offsets. 
  for (i = 1; i< VDimension; ++i)
    {
      accumulator += idx[i] * stride[i];
    }

  // Increment pointers.
  for (it = this->begin(); it < _end; ++it)
    {
      (*it) -= accumulator;
    }
  if (m_OutputBuffer)
    {
      m_OutputBuffer -= accumulator;
    }

  // Update loop counter values
  for (i=0; i<VDimension; ++i)
    {
      m_Loop[i]-= idx[i];
    }

  return *this;
}


} // namespace itk
