/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// It is important to use OrientedImages
#include "itkOrientedImage.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Project includes
#include "tubeARFFParser.h"

// Must include CLP before including tubeCLIHleperFunctions
#include "tubeWekaToImageCLP.h"

// No need for the DoIt nonsense (we arent' templated over input type)
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "WekaToImage",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef unsigned char                                 PixelType;
  const unsigned int                                    Dimension = 2;
  typedef itk::OrientedImage< PixelType,  Dimension >   ImageType;
  typedef itk::ImageFileWriter< ImageType >             WriterType;

  double progress = 0.1;
  progressReporter.Report( progress );

  tube::ARFFParser* parser = new tube::ARFFParser();
  parser->SetFilename( inputFile );
  timeCollector.Start("Load data");
  parser->Parse();
  timeCollector.Stop("Load data");

  std::list<float*> data = parser->GetARFFData();

  timeCollector.Start("Create image");
  ImageType::Pointer img = ImageType::New();
  ImageType::RegionType region;
  ImageType::SizeType size;
  ImageType::IndexType index;
  index[0] = parser->GetMinX();
  index[1] = parser->GetMinY();
  size[0] = (parser->GetMaxX() - parser->GetMinX())+1;
  size[1] = (parser->GetMaxY() - parser->GetMinY())+1;
  region.SetSize( size );
  region.SetIndex( index );
  img->SetRegions( region );
  img->Allocate();
  img->FillBuffer( 0 );

  for( std::list<float*>::const_iterator itr = data.begin();
       itr != data.end();
       ++itr )
    {
    ImageType::IndexType loc;
    loc[0] = (*itr)[0];
    loc[1] = (*itr)[1];
    img->SetPixel( loc, (*itr)[2] );
    }
  progress = 0.5;
  progressReporter.Report( progress );
  timeCollector.Stop("Create image");

  timeCollector.Start("Save data");
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( img );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Save data");

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  std::cout << std::endl;

  timeCollector.Report();

  delete parser;

  return EXIT_SUCCESS;
}