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
#include <QApplication>

// CTK includes
#include "ctkTest.h"

// MRML includes
#include "qSlicerLITTPlanV2Module.h"
#include "qSlicerLITTPlanV2ModuleWidget.h"
#include <vtkMRMLScene.h>
#include <vtkMRMLLinearTransformNode.h>

// VTK includes
#include <vtkMatrix4x4.h>
#include <vtkNew.h>

// ----------------------------------------------------------------------------
class qSlicerLITTPlanV2ModuleWidgetTester: public QObject
{
  Q_OBJECT

private slots:

  void testIdentity();
  void testInvert();
};

// ----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidgetTester::testIdentity()
{
  vtkNew<vtkMRMLScene> scene;
  vtkNew<vtkMRMLLinearTransformNode> transformNode;
  scene->AddNode(transformNode.GetPointer());

  qSlicerLITTPlanV2Module transformsModule;
  transformsModule.setMRMLScene(scene.GetPointer());
  transformsModule.logic();
  qSlicerLITTPlanV2ModuleWidget* transformsWidget =
    dynamic_cast<qSlicerLITTPlanV2ModuleWidget*>(transformsModule.widgetRepresentation());

  vtkMatrix4x4* matrix = transformNode->GetMatrixTransformToParent();
  matrix->SetElement(0,0, 10.);
  matrix->SetElement(1,0, 2.);
  transformsWidget->identity();
  QCOMPARE(matrix->GetElement(0,0), 1.);
  QCOMPARE(matrix->GetElement(1,0), 0.);
  //transformsWidget->show();
  //qApp->exec();
}

// ----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidgetTester::testInvert()
{
  vtkNew<vtkMRMLScene> scene;
  vtkNew<vtkMRMLLinearTransformNode> transformNode;
  scene->AddNode(transformNode.GetPointer());

  qSlicerLITTPlanV2Module transformsModule;
  transformsModule.setMRMLScene(scene.GetPointer());
  transformsModule.logic();
  qSlicerLITTPlanV2ModuleWidget* transformsWidget =
    dynamic_cast<qSlicerLITTPlanV2ModuleWidget*>(transformsModule.widgetRepresentation());

  vtkMatrix4x4* matrix = transformNode->GetMatrixTransformToParent();
  matrix->SetElement(0,0, 10.);
  matrix->SetElement(1,0, 2.);
  transformsWidget->invert();
  QCOMPARE(matrix->GetElement(0,0), 0.1);
  QCOMPARE(matrix->GetElement(1,0), -0.2);
  //transformsWidget->show();
  //qApp->exec();
}

// ----------------------------------------------------------------------------
CTK_TEST_MAIN(qSlicerLITTPlanV2ModuleWidgetTest)
#include "moc_qSlicerLITTPlanV2ModuleWidgetTest.cxx"
