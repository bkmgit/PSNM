// Adaptor for getting fortran simulation code into ParaView CoProcessor.
// Based on the PhastaAdaptor sample in the ParaView distribution.
// ParaView-3.14.1-Source/CoProcessing/Adaptors/FortranAdaptors/PhastaAdaptor/PhastaAdaptor.cxx


// Fortran specific header
// ParaView-3.14.1-Source/CoProcessing/Adaptors/FortranAdaptors/
#include "FortranAdaptorAPI.h" 

// CoProcessor specific headers
// ParaView-3.14.1-Source/CoProcessing/CoProcessor/
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"

#include "vtkDoubleArray.h"
#include "vtkPointData.h"

// ParaView-3.14.1-Source/VTK/Filtering/vtkImageData.h
#include "vtkImageData.h"

// These will be called from the Fortran "glue" code"
// Completely dependent on data layout, structured vs. unstructured, etc.
// since VTK/ParaView uses different internal layouts for each.

// Creates the data container for the CoProcessor.
extern "C" void createcpimagedata_(int* nx, int* ny, int* nz)
{
  if (!ParaViewCoProcessing::GetCoProcessorData()) {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
  }

  // The simulation grid is a 2-dimensional topologically and geometrically 
  // regular grid. In VTK/ParaView, this is considered an image data set.
  vtkImageData* Grid = vtkImageData::New();
  
  // assuming dimZ == 1 for now
  Grid->SetDimensions(*nx, *ny, *nz);
  
  // Setting the Origin and Spacing is also an option, but if the spacing is the same in each direction
  // probably not an issue. 

  // Name should be consistent between here, Fortran and Python client script.
  ParaViewCoProcessing::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(Grid);
}

// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// Might be an issue, VTK probably assumes row major, but
// omeg probably passed column major...
// by hand name mangling for fortran
extern "C" void addfield_(double* scalars)
{
  vtkCPInputDataDescription *idd = ParaViewCoProcessing::GetCoProcessorData()->GetInputDescriptionByName("input");

  // mvm: do I need this? Is this here because the client pipeline
  // grid might have been a different resolution?
  vtkImageData* Image = vtkImageData::SafeDownCast(idd->GetGrid());

  if (!Image) {
    vtkGenericWarningMacro("No adaptor grid to attach field data to.");
    return;
  }

  // vtkIdType NumberOfNodes = Image->GetNumberOfPoints();

  // assuming a single scalar, omega.
  // additional scalars or vectors would have similar setup
  // is omega being 2D going to be a problem...?

  // field name must match that in the fortran code.
  if (idd->IsFieldNeeded("omeg")) {
    vtkDoubleArray* omega = vtkDoubleArray::New();
    omega->SetName("omeg");
    omega->SetArray(scalars, Image->GetNumberOfPoints(), 1); // mvm 1 or 0? "save"?
    Image->GetPointData()->AddArray(omega);
    omega->Delete();
  }
}