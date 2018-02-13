/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSpeedFunctionToPathFilter.h,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkSpeedFunctionToPathFilter_h
#define itkSpeedFunctionToPathFilter_h

#include "itkArrivalFunctionToPathFilter.h"
#include "itkSpeedFunctionPathInformation.h"

namespace itk
{

/** \class SpeedFunctionToPathFilter
 * \brief Extracts a path from a speed function between a start point
 *       and end point, which also passes near the given way points.
 *
 * This filter propagates a number of Fast Marching fronts using the
 * given speed image to facilitate the path extraction.
 * For each end- and way-point, a Fast Marching arrival function is
 * computed. The minimal path is extracted for each segment and
 * concatenated to obtain the total path between the start- and
 * end-points, while passing near the given way-points.
 *
 * The user must provide the following:
 *    1. A real-valued speed function as the filter input
 *    2. At least one PathInfo object using AddPathInfo().
 * The speed function must be a real-valued (float or double) image
 * in the range [0,1]. If multiple PathInfo objects are added,
 * multiple paths are extracted and saved to separate filter outputs.
 *
 * A cost function optimizer may also be provided. If an optimizer
 * is not given, a RegularStepGradientDescentOptimizer is created
 * with default settings. Other suitable optimizers include
 * GradientDescentOptimizer and IterateNeighborhoodOptimizer.
 * See itkArrivalFunctionToPathFilter.h for more details.
 *
 * This filter is based on the methods described in:
 * [1] J. Sethian. Level Set Methods and Fast Marching Methods, chapter 20.
 *     Cambridge Press, 2nd edition, 1999.
 * [2] J. Andrews and J. Sethian. Fast marching methods for the continuous
 *     traveling salesman problem. Proceedings of the National Academy of
 *     Sciences (PNAS), 104(4):1118 1123, 2007.
 *
 * \author Dan Mueller, Queensland University of Technology, dan.muel[at]gmail.com
 *
 * \sa ArrivalFunctionToPathFilter
 * \ingroup ImageToPathFilters
 *
 * \ingroup MinimalPathExtraction
 */
template <class TInputImage,
          class TOutputPath = PolyLineParametricPath<TInputImage::ImageDimension> >
class ITK_EXPORT SpeedFunctionToPathFilter :
    public ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
{
public:
  /** Standard class type alias. */
  using Self = SpeedFunctionToPathFilter;
  using Superclass = ArrivalFunctionToPathFilter<TInputImage,TOutputPath>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(SpeedFunctionToPathFilter,ArrivalFunctionToPathFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** ImageDimension constants */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;

  /** Some image type alias. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using InputImageRegionType = typename InputImageType::RegionType;
  using InputImagePixelType = typename InputImageType::PixelType;

  /** Some path type alias. */
  using OutputPathType = TOutputPath;
  using OutputPathPointer = typename OutputPathType::Pointer;
  using OutputPathConstPointer = typename OutputPathType::ConstPointer;

  /** Some convenient type alias. */
  using ContinuousIndexType = typename Superclass::ContinuousIndexType;
  using IndexType = typename Superclass::IndexType;
  using PointType = typename Superclass::PointType;
  using CostFunctionType = typename Superclass::CostFunctionType;
  using OptimizerType = typename Superclass::OptimizerType;

  /** Path information type alias. */
  using PathInformationType = SpeedFunctionPathInformation<PointType>;
  using PointsContainerType = typename PathInformationType::PointsContainerType;

  /** Override superclass behaviour.
   *  Warning: SetPathEndPoint() is not valid for this filter.
   *  This method is provided by the superclass, however it is not
   *  used by this subclass. Use AddPathInfo() instead.*/
  void SetPathEndPoint( const PointType& ) override
  {
    itkWarningMacro("SetPathEndPoint() is not valid for this filter. Use AddPathInfo() instead.");
  }

  /** Override superclass behaviour.
   *  Warning: AddPathEndPoint() is not valid for this filter.
   *  This method is provided by the superclass, however it is not
   *  used by this subclass. Use AddPathInfo() instead.*/
  void AddPathEndPoint( const PointType& ) override
  {
    itkWarningMacro("AddPathEndPoint() is not valid for this filter. Use AddPathInfo() instead.");
  }

  /** Override superclass behaviour.
   *  Warning: ClearPathEndPoints() is not valid for this filter.
   *  This method is provided by the superclass, however it is not
   *  used by this subclass. Use ClearPathInfo() instead.*/
  void ClearPathEndPoints() override
  {
    itkWarningMacro("ClearPathEndPoints() is not valid for this filter. Use ClearPathInfo() instead.");
  }

  /** Add a path information object to process.
   *  At least one PathInfo object must be added before processing. */
  void AddPathInformation(PathInformationType * info )
  {
    m_Information.push_back( info );
  }

  /** Clear the list of path information objects. */
  void ClearPathInformation()
  {
    m_Information.clear( );
  }

  /** Handle optimizer iteration events. */
  void Execute( const itk::Object * object, const itk::EventObject & event ) override;

  /** access the arrival image for debugging purposes */
  itkGetConstMacro( CurrentArrivalFunction, InputImagePointer );

protected:
  SpeedFunctionToPathFilter( );
  ~SpeedFunctionToPathFilter( ) override;
  void PrintSelf( std::ostream& os, Indent indent ) const override;

  /** Implemention of algorithm. */
  void GenerateData( void ) override;

  /** Get the number of paths which the user has instructed to extracted. */
  unsigned int GetNumberOfPathsToExtract( ) const override;

  /** Compute the arrival function from which to extract the path. */
  InputImageType * ComputeArrivalFunction( ) override;

  /** Override handling of optimizer iteration events to accomodate way points. */
  const PointsContainerType & GetNextEndPoint( ) override;

  std::vector< typename PathInformationType::Pointer > m_Information;
  InputImagePointer                                    m_CurrentArrivalFunction;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(SpeedFunctionToPathFilter);

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpeedFunctionToPathFilter.hxx"
#endif

#endif
