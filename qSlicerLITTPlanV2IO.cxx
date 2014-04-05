/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Julien Finet, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// Qt includes

// SlicerQt includes
#include "qSlicerLITTPlanV2IO.h"

// Logic includes
#include "vtkSlicerTransformLogic.h"

// MRML includes
#include <vtkMRMLTransformNode.h>

// VTK includes
#include <vtkSmartPointer.h>

//-----------------------------------------------------------------------------
class qSlicerLITTPlanV2IOPrivate
{
public:
  vtkSmartPointer<vtkSlicerTransformLogic> TransformLogic;
};

//-----------------------------------------------------------------------------
qSlicerLITTPlanV2IO::qSlicerLITTPlanV2IO(
  vtkSlicerTransformLogic* _transformLogic, QObject* _parent)
  : qSlicerIO(_parent)
  , d_ptr(new qSlicerLITTPlanV2IOPrivate)
{
  this->setTransformLogic(_transformLogic);
}


//-----------------------------------------------------------------------------
qSlicerLITTPlanV2IO::~qSlicerLITTPlanV2IO()
{
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2IO::setTransformLogic(vtkSlicerTransformLogic* newTransformLogic)
{
  Q_D(qSlicerLITTPlanV2IO);
  d->TransformLogic = newTransformLogic;
}

//-----------------------------------------------------------------------------
vtkSlicerTransformLogic* qSlicerLITTPlanV2IO::transformLogic()const
{
  Q_D(const qSlicerLITTPlanV2IO);
  return d->TransformLogic;
}

//-----------------------------------------------------------------------------
QString qSlicerLITTPlanV2IO::description()const
{
  return "Transform";
}

//-----------------------------------------------------------------------------
qSlicerIO::IOFileType qSlicerLITTPlanV2IO::fileType()const
{
  return QString("TransformFile");
}

//-----------------------------------------------------------------------------
QStringList qSlicerLITTPlanV2IO::extensions()const
{
  return QStringList() << "Transform (*.tfm *.mat *.txt)";
}

//-----------------------------------------------------------------------------
bool qSlicerLITTPlanV2IO::load(const IOProperties& properties)
{
  Q_D(qSlicerLITTPlanV2IO);
  Q_ASSERT(properties.contains("fileName"));
  QString fileName = properties["fileName"].toString();

  if (d->TransformLogic.GetPointer() == 0)
    {
    return false;
    }
  vtkMRMLTransformNode* node = d->TransformLogic->AddTransform(
    fileName.toLatin1(), this->mrmlScene());
  if (node)
    {
   // this->setLoadedNodes(QStringList(QString(node->GetID())));
    }
  else
    {
    //this->setLoadedNodes(QStringList());
    }
  return node != 0;
}

// TODO: add the save() method. Use vtkSlicerTransformLogic::SaveTransform()
