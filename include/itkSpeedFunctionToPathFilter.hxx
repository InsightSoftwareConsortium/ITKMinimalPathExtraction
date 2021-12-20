/*=========================================================================
 *
 *  Copyright NumFOCUS
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
#ifndef itkSpeedFunctionToPathFilter_hxx
#define itkSpeedFunctionToPathFilter_hxx

#include "itkMath.h"
#include "itkFastMarchingUpwindGradientImageFilter.h"
#include <itkConstNeighborhoodIterator.h>
#include <itkConstantBoundaryCondition.h>

namespace itk
{

template <typename TInputImage, typename TOutputPath>
SpeedFunctionToPathFilter<TInputImage, TOutputPath>::SpeedFunctionToPathFilter()
  : m_CurrentArrivalFunction(nullptr)
{
  m_TargetRadius = 2;
  m_AutoTerminate = true;
  m_AutoTerminateFactor = 0.5;
}


template <typename TInputImage, typename TOutputPath>
SpeedFunctionToPathFilter<TInputImage, TOutputPath>::~SpeedFunctionToPathFilter() = default;


template <typename TInputImage, typename TOutputPath>
unsigned int
SpeedFunctionToPathFilter<TInputImage, TOutputPath>::GetNumberOfPathsToExtract() const
{
  return m_Information.size();
}


template <typename TInputImage, typename TOutputPath>
const typename SpeedFunctionToPathFilter<TInputImage, TOutputPath>::PointsContainerType &
SpeedFunctionToPathFilter<TInputImage, TOutputPath>::GetNextEndPoint()
{
  return m_Information[Superclass::m_CurrentOutput]->GetEndPoint();
}

/**
* take a collection of indexes and produce the equivalent of a
* dilated version.
*/

template<typename TInputImage, typename TOutputPath>
typename SpeedFunctionToPathFilter<TInputImage,TOutputPath>::IndexTypeSet
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::GetNeighbors(IndexTypeVec idxs)
{
  InputImagePointer speed =
    const_cast< InputImageType * >( this->GetInput() );

  using BoundaryConditionType = ConstantBoundaryCondition<InputImageType>;

  IndexTypeSet UniqueIndexes;
  typename InputImageType::SizeType radius;
  radius.Fill(m_TargetRadius);
  ConstNeighborhoodIterator<InputImageType, BoundaryConditionType>
    niterator(radius, speed, speed->GetLargestPossibleRegion());

  BoundaryConditionType bc;
  bc.SetConstant(0);
  niterator.SetBoundaryCondition(bc);
  niterator.NeedToUseBoundaryConditionOn();

  for (auto it = idxs.begin(); it != idxs.end(); it++ )
    {
      niterator.SetLocation(*it);
      if ( niterator.GetCenterPixel() > 0 )
	{
	  // Visit the entire neighborhood (including center) and
	  // add any pixel that has a nonzero speed function value
	  for (auto NB = 0; NB < niterator.Size(); NB++)
	    {
	      if ( niterator.GetPixel(NB) > 0 )
		{
		  UniqueIndexes.insert(niterator.GetIndex(NB));
		}
	    }
	}
    }

  return(UniqueIndexes);
}


template <typename TInputImage, typename TOutputPath>
typename SpeedFunctionToPathFilter<TInputImage,TOutputPath>::InputImagePixelType
SpeedFunctionToPathFilter<TInputImage,TOutputPath>::GetTrialGradient(IndexTypeVec idxs)
{
  InputImagePointer arrival =  m_CurrentArrivalFunction;

  using BoundaryConditionType = ConstantBoundaryCondition<InputImageType>;

  IndexTypeSet UniqueIndexes;
  typename InputImageType::SizeType radius;
  radius.Fill(1);
  ConstNeighborhoodIterator<InputImageType, BoundaryConditionType>
    niterator(radius, arrival, arrival->GetLargestPossibleRegion());

  BoundaryConditionType bc;
  bc.SetConstant(itk::NumericTraits<InputImagePixelType>::max());
  niterator.SetBoundaryCondition(bc);
  niterator.NeedToUseBoundaryConditionOn();

  // looking for the smallest nonzero difference
  InputImagePixelType mindiff(itk::NumericTraits<InputImagePixelType>::max());

  for (auto it = idxs.begin(); it != idxs.end(); it++ )
    {
      niterator.SetLocation(*it);
      InputImagePixelType CP = niterator.GetCenterPixel();
      // Visit the entire neighborhood (including center) and
      // add any pixel that has a nonzero arrival function value
      for (auto NB = 0; NB < niterator.Size(); NB++)
	{
	  // CP values should always be zero
	  InputImagePixelType NPD = niterator.GetPixel(NB) - CP;
	  if (NPD  > 0 )
	    {
	      mindiff=std::min(mindiff, NPD);
	    }
	}
    }

  return(mindiff);
}


template <typename TInputImage, typename TOutputPath>
typename SpeedFunctionToPathFilter<TInputImage, TOutputPath>::InputImageType *
SpeedFunctionToPathFilter<TInputImage, TOutputPath>::ComputeArrivalFunction()
{
  // Get the speed image
  InputImagePointer speed = const_cast<InputImageType *>(this->GetInput());

  // Set the fast marching method for computing the arrival function
  using FastMarchingType = FastMarchingUpwindGradientImageFilter<TInputImage, TInputImage>;

  using NodeContainer = typename FastMarchingType::NodeContainer;
  using NodeType = typename FastMarchingType::NodeType;
  typename FastMarchingType::Pointer marching = FastMarchingType::New();
  marching->SetInput(speed);
  marching->SetGenerateGradientImage(false);
  marching->SetTargetReachedModeToAllTargets();
  // Add next and previous front sources as target points to
  // limit the front propagation to just the required zones
  PointsContainerType PrevFront = m_Information[Superclass::m_CurrentOutput]->PeekPreviousFront();
  PointsContainerType NextFront = m_Information[Superclass::m_CurrentOutput]->PeekNextFront();
  IndexTypeVec PrevIndexVec(0);
  IndexTypeVec NextIndexVec(0);

  typename NodeContainer::Pointer targets = NodeContainer::New();
  targets->Initialize();

  for (auto it = PrevFront.begin(); it != PrevFront.end(); it++)
  {
    IndexType indexTargetPrevious;
    NodeType  nodeTargetPrevious;
    speed->TransformPhysicalPointToIndex(*it, indexTargetPrevious);
    PrevIndexVec.push_back(indexTargetPrevious);
  }

  for (auto it = NextFront.begin(); it != NextFront.end(); it++)
  {
    IndexType indexTargetNext;
    NodeType  nodeTargetNext;
    speed->TransformPhysicalPointToIndex(*it, indexTargetNext);
    NextIndexVec.push_back(indexTargetNext);
  }
  IndexTypeVec AllTargets( PrevIndexVec );
  AllTargets.insert(AllTargets.end(), NextIndexVec.begin(), NextIndexVec.end());

  IndexTypeSet UniqueTargets = GetNeighbors( AllTargets );
  // Add neighbours of all the targets points to ensure that the
  // gradients in the neighborhood of each potential destination point
  // is smooth.
  for (auto it = UniqueTargets.begin(); it != UniqueTargets.end(); it++)
    {
      NodeType nodeTarget;
      nodeTarget.SetValue( 0.0 );
      nodeTarget.SetIndex( *it );
      targets->InsertElement( targets->Size(), nodeTarget );
    }
  marching->SetTargetPoints( targets );

  // Get the next Front source point and add as trial point
  typename NodeContainer::Pointer trial = NodeContainer::New();
  trial->Initialize();
  PointsContainerType CurrentFront =
    m_Information[Superclass::m_CurrentOutput]->PeekCurrentFront(); // FrontAndAdvance();
  IndexTypeVec CurrentIndexVec(0);

  for (auto it = CurrentFront.begin(); it != CurrentFront.end(); it++)
  {
    IndexType indexTrial;
    NodeType  nodeTrial;
    speed->TransformPhysicalPointToIndex(*it, indexTrial);
    nodeTrial.SetValue(0.0);
    nodeTrial.SetIndex(indexTrial);
    trial->InsertElement(trial->Size(), nodeTrial);
    CurrentIndexVec.push_back(indexTrial);
  }
  marching->SetTrialPoints(trial);

  // Update the method and set the arrival function
  marching->UpdateLargestPossibleRegion();
  m_CurrentArrivalFunction = marching->GetOutput();
  m_CurrentArrivalFunction->DisconnectPipeline();

  // Only the index with the minimum arrival time should stay in the "Previous" point set
  // This will be used to initialise the optimizer
  if (PrevFront.size() > 1)
  {
    InputImagePixelType MinTime = itk::NumericTraits<InputImagePixelType>::max();
    unsigned            MinPos(0);
    for (unsigned idx = 0; idx < PrevIndexVec.size(); ++idx)
    {
      InputImagePixelType V = m_CurrentArrivalFunction->GetPixel(PrevIndexVec[idx]);
      if (V < MinTime)
      {
        MinPos = idx;
        MinTime = V;
      }
    }
    m_Information[Superclass::m_CurrentOutput]->SetPrevious(PrevFront[MinPos]);
  }

  if (m_AutoTerminate)
    {
      // Examine the neighbours of the trial points to determine the minimum neighbour
      // difference, for the purpose of estimating a good termination value.
      InputImagePixelType MD = GetTrialGradient( CurrentIndexVec );
      std::cout << "Min diff for termination = " << MD << std::endl;
      this->SetTerminationValue( MD * this->GetAutoTerminateFactor() );
    }
  // Make the arrival function flat inside the seeds, otherwise the
  // optimizer will cross over them. This only matters if the seeds are extended
  // Not sure that this is needed any more. Expect it was a side effect of only
  // adding a single point to the trial set.
  if (CurrentIndexVec.size() > 1)
  {
    for (auto vi = CurrentIndexVec.begin(); vi != CurrentIndexVec.end(); vi++)
    {
      m_CurrentArrivalFunction->SetPixel(*vi, 0);
    }
  }

  m_Information[Superclass::m_CurrentOutput]->Advance();
  return m_CurrentArrivalFunction;
}


template <typename TInputImage, typename TOutputPath>
void
SpeedFunctionToPathFilter<TInputImage, TOutputPath>::GenerateData()
{
  // Get the speed function
  InputImagePointer speed = const_cast<InputImageType *>(this->GetInput());
  if (speed.IsNull())
  {
    itkExceptionMacro("Speed function image must be provided");
  }

  // Ensure the user has added at least one path info object
  if (m_Information.empty())
  {
    itkExceptionMacro("No PathInfo objects: at least one must be added.");
  }

  // Extract the path
  Superclass::GenerateData();
}


template <typename TInputImage, typename TOutputPath>
void
SpeedFunctionToPathFilter<TInputImage, TOutputPath>::Execute(const Object *      object,
                                                             const EventObject & itkNotUsed(event))
{
  // Cast object to optmizer
  typename OptimizerType::Pointer optimizer = (OptimizerType *)dynamic_cast<const OptimizerType *>(object);
  if (optimizer.IsNull())
    return;

  // Get current position and value
  typename OptimizerType::ParametersType currentParameters = optimizer->GetCurrentPosition();
  unsigned int                           lenParameters = currentParameters.GetSize();
  if (lenParameters != InputImageDimension)
    return;
  typename OptimizerType::MeasureType currentValue = optimizer->GetValue(currentParameters);

  // Convert parameters to point
  bool         valid = false;
  unsigned int numparams = optimizer->GetCurrentPosition().GetSize();
  PointType    point;
  point.Fill(0.0);
  for (unsigned int i = 0; i < numparams; i++)
  {
    point[i] = optimizer->GetCurrentPosition()[i];
    valid = true;
  }
  if (!valid)
    return;


  // Check if we have reached the termination value
  if (currentValue < Superclass::m_TerminationValue && m_Information[Superclass::m_CurrentOutput]->HasNextFront())
  {
    // We have terminated the current path segment,
    // but there are more fronts to propagate

    // TODO: The path has not actually reached the path point.
    //       Change the next front point to be the current point.
    // only the arrival point reached by the optimizer should be included in
    // the extended point set
    if (m_Information[Superclass::m_CurrentOutput]->PeekPreviousFront().size() > 1)
    {
      m_Information[Superclass::m_CurrentOutput]->SetPrevious(point);
    }
    // Update the arrival function and re-initialise the cost function
    Superclass::m_CostFunction->SetImage(this->ComputeArrivalFunction());
    Superclass::m_CostFunction->Initialize();
  }
  else if (currentValue >= Superclass::m_TerminationValue)
  {
    // Convert point to continuous index
    InputImagePointer   input = const_cast<InputImageType *>(this->GetInput());
    ContinuousIndexType cindex;
    input->TransformPhysicalPointToContinuousIndex(point, cindex);

    // Add point as vertex in path
    OutputPathPointer output = this->GetOutput(Superclass::m_CurrentOutput);
    output->AddVertex(cindex);
  }
}


template <typename TInputImage, typename TOutputPath>
void
SpeedFunctionToPathFilter<TInputImage, TOutputPath>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


} // end namespace itk

#endif
