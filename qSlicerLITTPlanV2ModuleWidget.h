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

  This file was created and modified by Erol Yeniaras using the template originally 
  developed by Jean-Christophe Fillion-Robin.

==============================================================================*/
#ifndef __qSlicerLITTPlanV2ModuleWidget_h
#define __qSlicerLITTPlanV2ModuleWidget_h

// CTK includes
#include <ctkPimpl.h>

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

// LITTPlanV2 includes
#include "qSlicerLITTPlanV2ModuleExport.h"
#include <QFileSystemWatcher>

#include <vtkCallbackCommand.h>
#include <vtkSmartPointer.h>

class vtkMatrix4x4;
class vtkMRMLNode;
class qSlicerLITTPlanV2ModuleWidgetPrivate;
class vtkMRMLAnnotationFiducialNode;
class vtkMRMLFiducialListNode;
class vtkMRMLAnnotationHierarchyNode;

class Q_SLICER_QTMODULES_LITTPLANV2_EXPORT qSlicerLITTPlanV2ModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT
  QVTK_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerLITTPlanV2ModuleWidget(QWidget *parent=0);
  virtual ~qSlicerLITTPlanV2ModuleWidget();

  /// Reimplemented for internal reasons
  void setMRMLScene(vtkMRMLScene* scene);
  
  vtkSmartPointer<vtkCallbackCommand> CallBack1; // Observe fiducial 2
  vtkSmartPointer<vtkCallbackCommand> CallBack2; // Observe fiducial 1

  vtkSmartPointer<vtkMRMLAnnotationFiducialNode> fnode0; // Start point of the applicator path
  vtkSmartPointer<vtkMRMLAnnotationFiducialNode> fnode1; // End point of the applicator path

  QFileSystemWatcher *femWatcher;

  bool observersAreActive;
  bool fiducialsAreValid;
  bool applicatorVisible;
  bool pathVisible;  
  double wmtParameters[5]; // White matter parameters
  double gmtParameters[5]; // Grey matter parameters
  double csfParameters[5]; // CSF tissue parameters
  double tumParameters[5]; // Tumor tissue parameters
  double vesParameters[5]; // Vessel parameters
  double edeParameters[5]; // Vessel parameters
  double othParameters[5]; // Other tissue parameters
  double power; // A global variable that is same for every tissue
  double thermalDose; // A global variable that is same for every tissue == This is Isoterm now!!! (Updated 08/14-2012)
  int nelemtip; // A global variable
  double powerCoefficient; // A global variable

  double radius; // A global variable
  int tissueIndex; // Selected tissue for setting the parameters and appying the treatment
  double startPoint[3]; // Fiducial points coordinates for guiding the intervention
  double targetPoint[3]; // They are entrance and target points respectively
  double prevStartPoint[3]; // Previous coordinates of fiducial points
  double prevTargetPoint[3]; // If they were created previously these points are full otherwise null
  
public slots:
   
	void onFiducalNodeChanged(vtkMRMLNode* node);
    void LoadLabelMap(); // Load the initial label map
	void CreatePath(); // Create the applicator path
	void CreatePathOld(); // Create the applicator path
	void CreateSphereAtTarget(); // Create the applicator path
    void LoadMesh();
	void HidePath();
	void LoadApplicator(); // Place applicator on the path
	void sleep(unsigned int mseconds);
  
  void LoadParametersToGUI(); // Load parameters to textboxes to show them on the GUI
  void LoadParametersFromGUI(); // Load parameters from GUI to the arrays
  void SetParametersToDefault(); // Set the parameters to defaults
  
  void onBtnDefaultClicked();
  void onTissueTypeChanged();
  //void GetFEM();
  //void onBtnPathClicked();
  void onBtnApplicatorClicked();
  void onBtnBurnClicked();
  static void OnFiducial1Moved(vtkObject* vtk_obj, unsigned long event, void* client_data, void* call_data);
  static void OnFiducial2Moved(vtkObject* vtk_obj, unsigned long event, void* client_data, void* call_data);

  /// Set the matrix to identity, the sliders are reset to the position 0
  void identity();

  /// Invert the matrix. The sliders are reset to the position 0.
  void invert();

protected:
  virtual void setup();

protected slots:
  
  void onCoordinateReferenceButtonPressed(int id);
  void onNodeSelected(vtkMRMLNode* node);
  void onTranslationRangeChanged(double newMin, double newMax);

  void transformSelectedNodes();
  void untransformSelectedNodes();
  /// 
  /// Triggered upon MRML transform node updates
  void onMRMLTransformNodeModified(vtkObject* caller);

protected:
  /// 
  /// Fill the 'minmax' array with the min/max translation value of the matrix.
  /// Parameter expand allows to specify (using a value between 0 and 1)
  /// which percentage of the found min/max value should be substracted/added
  /// to the min/max value found.
  void extractMinMaxTranslationValue(vtkMatrix4x4 * mat, double& min, double& max);

  /// 
  /// Convenient method to return the coordinate system currently selected
  int coordinateReference()const;
  bool EqualD(double a, double b);
  bool EqualP(double p1[], double p2[]);

protected:
  QScopedPointer<qSlicerLITTPlanV2ModuleWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerLITTPlanV2ModuleWidget);
  Q_DISABLE_COPY(qSlicerLITTPlanV2ModuleWidget);
};

#endif
