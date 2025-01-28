// General includes
#include <string>
#include <iostream>

// ITK includes
#include "itkNumericTraits.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPolyLineParametricPath.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkArrivalFunctionToPathFilter.h"
#include "itkSpeedFunctionToPathFilter.h"
#include "itkPathIterator.h"
#include "itkGradientDescentOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkIterateNeighborhoodOptimizer.h"

int
example_gradientdescent(int argc, char * argv[])
{
  // Typedefs
  constexpr unsigned int Dimension = 2;
  using PixelType = float;
  using OutputPixelType = unsigned char;
  using ImageType = itk::Image<PixelType, Dimension>;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;
  using ReaderType = itk::ImageFileReader<ImageType>;
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  using PathType = itk::PolyLineParametricPath<Dimension>;
  using PathFilterType = itk::SpeedFunctionToPathFilter<ImageType, PathType>;
  using CoordinateType = PathFilterType::CostFunctionType::CoordinateType;
  using PathIteratorType = itk::PathIterator<OutputImageType, PathType>;

  // Get filename arguments
  unsigned int argi = 1;
  const char * outputFileName = argv[argi++];
  const char * speedFileName = argv[argi++];

  // Read speed function
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(speedFileName);
  reader->Update();
  ImageType::Pointer speed = reader->GetOutput();
  speed->DisconnectPipeline();

  // Create interpolator
  using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, CoordinateType>;
  InterpolatorType::Pointer interp = InterpolatorType::New();

  // Create cost function
  PathFilterType::CostFunctionType::Pointer cost = PathFilterType::CostFunctionType::New();
  cost->SetInterpolator(interp);

  // Create optimizer
  using OptimizerType = itk::GradientDescentOptimizer;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetNumberOfIterations(1000);

  // Create path filter
  PathFilterType::Pointer pathFilter = PathFilterType::New();
  pathFilter->SetInput(speed);
  pathFilter->SetCostFunction(cost);
  pathFilter->SetOptimizer(optimizer);
  pathFilter->SetTerminationValue(2.0);

  // Setup path points
  PathFilterType::PointType start, end, way1;
  start[0] = 10;
  start[1] = 100;
  end[0] = 100;
  end[1] = 10;
  way1[0] = 10;
  way1[1] = 10;

  // Add path information
  using PathInformationType = PathFilterType::PathInformationType;
  PathInformationType::Pointer info = PathInformationType::New();
  info->SetStartPoint(start);
  info->SetEndPoint(end);
  info->AddWayPoint(way1);
  pathFilter->AddPathInformation(info);

  // Compute the path
  pathFilter->Update();

  // Allocate output image
  OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions(speed->GetLargestPossibleRegion());
  output->SetSpacing(speed->GetSpacing());
  output->SetOrigin(speed->GetOrigin());
  output->Allocate();
  output->FillBuffer(itk::NumericTraits<OutputPixelType>::Zero);

  // Rasterize path
  for (unsigned int i = 0; i < pathFilter->GetNumberOfOutputs(); i++)
  {
    // Get the path
    PathType::Pointer path = pathFilter->GetOutput(i);

    // Check path is valid
    if (path->GetVertexList()->Size() == 0)
    {
      std::cout << "WARNING: Path " << (i + 1) << " contains no points!" << std::endl;
      continue;
    }

    // Iterate path and convert to image
    PathIteratorType it(output, path);
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      it.Set(itk::NumericTraits<OutputPixelType>::max());
    }
  }

  // Write output
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFileName);
  writer->SetInput(output);
  writer->Update();

  // Return
  return EXIT_SUCCESS;
}

int
example_regularstepgradientdescent(int argc, char * argv[])
{
  // Typedefs
  constexpr unsigned int Dimension = 2;
  using PixelType = float;
  using OutputPixelType = unsigned char;
  using ImageType = itk::Image<PixelType, Dimension>;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;
  using ReaderType = itk::ImageFileReader<ImageType>;
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  using PathType = itk::PolyLineParametricPath<Dimension>;
  using PathFilterType = itk::SpeedFunctionToPathFilter<ImageType, PathType>;
  using CoordinateType = PathFilterType::CostFunctionType::CoordinateType;
  using PathIteratorType = itk::PathIterator<OutputImageType, PathType>;

  // Get filename arguments
  unsigned int argi = 1;
  const char * outputFileName = argv[argi++];
  const char * speedFileName = argv[argi++];

  // Read speed function
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(speedFileName);
  reader->Update();
  ImageType::Pointer speed = reader->GetOutput();
  speed->DisconnectPipeline();

  // Create interpolator
  using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, CoordinateType>;
  InterpolatorType::Pointer interp = InterpolatorType::New();

  // Create cost function
  PathFilterType::CostFunctionType::Pointer cost = PathFilterType::CostFunctionType::New();
  cost->SetInterpolator(interp);

  // Create optimizer
  using OptimizerType = itk::RegularStepGradientDescentOptimizer;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetNumberOfIterations(1000);
  optimizer->SetMaximumStepLength(0.5);
  optimizer->SetMinimumStepLength(0.1);
  optimizer->SetRelaxationFactor(0.5);

  // Create path filter
  PathFilterType::Pointer pathFilter = PathFilterType::New();
  pathFilter->SetInput(speed);
  pathFilter->SetCostFunction(cost);
  pathFilter->SetOptimizer(optimizer);
  pathFilter->SetTerminationValue(2.0);

  // Setup path points
  PathFilterType::PointType start, end, way1;
  start[0] = 10;
  start[1] = 100;
  end[0] = 100;
  end[1] = 10;
  way1[0] = 10;
  way1[1] = 10;

  // Add path information
  using PathInformationType = PathFilterType::PathInformationType;
  PathInformationType::Pointer info = PathInformationType::New();
  info->SetStartPoint(start);
  info->SetEndPoint(end);
  info->AddWayPoint(way1);
  pathFilter->AddPathInformation(info);

  // Compute the path
  pathFilter->Update();

  // Allocate output image
  OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions(speed->GetLargestPossibleRegion());
  output->SetSpacing(speed->GetSpacing());
  output->SetOrigin(speed->GetOrigin());
  output->Allocate();
  output->FillBuffer(itk::NumericTraits<OutputPixelType>::Zero);

  // Rasterize path
  for (unsigned int i = 0; i < pathFilter->GetNumberOfOutputs(); i++)
  {
    // Get the path
    PathType::Pointer path = pathFilter->GetOutput(i);

    // Check path is valid
    if (path->GetVertexList()->Size() == 0)
    {
      std::cout << "WARNING: Path " << (i + 1) << " contains no points!" << std::endl;
      continue;
    }

    // Iterate path and convert to image
    PathIteratorType it(output, path);
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      it.Set(itk::NumericTraits<OutputPixelType>::max());
    }
  }

  // Write output
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFileName);
  writer->SetInput(output);
  writer->Update();

  // Return
  return EXIT_SUCCESS;
}

int
main(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: MinimalPathExamples <SpeedImage> <OutputImage> [UseRegularStepGradientDescentOptimizer]"
              << std::endl;
    return EXIT_FAILURE;
  }

  int useRegularStepGradientDescentOptimizer = 0;
  if (argc > 3)
  {
    useRegularStepGradientDescentOptimizer = std::atoi(argv[3]);
  }

  if (useRegularStepGradientDescentOptimizer)
  {
    return example_regularstepgradientdescent(argc, argv);
  }
  return example_gradientdescent(argc, argv);
}
