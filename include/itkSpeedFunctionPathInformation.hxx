/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSpeedFunctionPathInformation.hxx,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef itkSpeedFunctionPathInformation_hxx
#define itkSpeedFunctionPathInformation_hxx

#include "itkSpeedFunctionPathInformation.h"

namespace itk
{

template <typename TPoint>
SpeedFunctionPathInformation<TPoint>
::SpeedFunctionPathInformation()
{
  this->ClearInfo();
}


template <typename TPoint>
SpeedFunctionPathInformation<TPoint>
::~SpeedFunctionPathInformation()
{
}


template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::ClearInfo()
{
  m_Information.clear();
  m_Information.resize(2);
  m_Front = 1;
}


template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetStartPoint( const PointType & start )
{
  m_Information[1] = PtoPVec(start);
}


template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetEndPoint( const PointType & end )
{
  m_Information[0] = PtoPVec( end );
}


template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::AddWayPoint( const PointType & way )
{
  m_Information.push_back( PtoPVec(way) );
  m_Front++;
}
  
template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetCurrent( const PointType & newcurrent )
{
  m_Information[m_Front] = PtoPVec(newcurrent);
}

template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::Advance( void )
{
  m_Front--;
}

template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetStartPoint( const PointsContainerType & start )
{
  m_Information[1] = start;
}


template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetEndPoint( const PointsContainerType & end )
{
  m_Information[0] = end;
}


template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::AddWayPoint( const PointsContainerType & way )
{
  m_Information.push_back( way );
  m_Front++;
}

template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetCurrent( const PointsContainerType & newcurrent )
{
  m_Information[m_Front] = newcurrent;
}
  
template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetPrevious( const PointsContainerType & newprevious )
{
   SizeValueType F;
  if ( m_Front == m_Information.size()-1 )
  {
    F=0;
  }
  else
  {
    F = m_Front+1;
  }
  m_Information[F] = newprevious;
}
template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetNext( const PointsContainerType & newnext )
{
  SizeValueType F;
  if (m_Front <= 1)
    {
      F=1;
    }
  else
    {
      F = m_Front - 1;
    }
  m_Information[F] = newnext;
}

template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetPrevious( const PointType & newprevious )
{
  SetPrevious(PtoPVec(newprevious));
}

template <typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::SetNext( const PointType & newnext )
{
  SetNext(PtoPVec(newnext));
}
  

template <typename TPoint>
unsigned int SpeedFunctionPathInformation<TPoint>
::GetNumberOfPoints( ) const
{
  return m_Information.size();
}


template <typename TPoint>
const typename SpeedFunctionPathInformation<TPoint>::PointsContainerType &
SpeedFunctionPathInformation<TPoint>
::GetStartPoint( ) const
{
  return m_Information[1];
}


template <typename TPoint>
const typename SpeedFunctionPathInformation<TPoint>::PointsContainerType &
SpeedFunctionPathInformation<TPoint>
::GetEndPoint( ) const
{
  return m_Information[0];
}


template <typename TPoint>
const typename SpeedFunctionPathInformation<TPoint>::PointsContainerType &
SpeedFunctionPathInformation<TPoint>
::GetWayPoint( SizeValueType i ) const
{
  return m_Information[2+i];
}


template <typename TPoint>
bool SpeedFunctionPathInformation<TPoint>
::HasNextFront( ) const
{
  return m_Front >= 1;

}


template <typename TPoint>
const typename SpeedFunctionPathInformation<TPoint>::PointsContainerType &
SpeedFunctionPathInformation<TPoint>
::GetCurrentFrontAndAdvance( )
{
  return m_Information[m_Front--];
}


template <typename TPoint>
const typename SpeedFunctionPathInformation<TPoint>::PointsContainerType &
SpeedFunctionPathInformation<TPoint>
::PeekCurrentFront( ) const
{
  return m_Information[m_Front];
}


template <typename TPoint>
const typename SpeedFunctionPathInformation<TPoint>::PointsContainerType &
SpeedFunctionPathInformation<TPoint>
::PeekNextFront( ) const
{
  if ( m_Front <= 1 )
  {
  return m_Information[1];
  }
  else
  {
  return m_Information[m_Front-1];
  }
}


template <typename TPoint>
const typename SpeedFunctionPathInformation<TPoint>::PointsContainerType &
SpeedFunctionPathInformation<TPoint>
::PeekPreviousFront( ) const
{
  if ( m_Front == m_Information.size()-1 )
  {
  return m_Information[0];
  }
  else
  {
  return m_Information[m_Front+1];
  }
}


template<typename TPoint>
void SpeedFunctionPathInformation<TPoint>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
