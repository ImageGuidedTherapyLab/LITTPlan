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

  This file was created and modified by Erol Yeniaras using the template (Transforms)
  developed by Jean-Christophe Fillion-Robin.

==============================================================================*/
#include "ctkUtils.h"
#include <fstream>
#include <time.h>

// Qt includes
#include <QFileDialog>
#include <QFileSystemWatcher>

// SlicerQt includes
#include "qSlicerLITTPlanV2ModuleWidget.h"
#include "ui_qSlicerLITTPlanV2Module.h"
#include "qSlicerApplication.h"
#include "qSlicerIOManager.h"

// vtkSlicerLogic includes
#include "vtkSlicerTransformLogic.h"

// MRMLWidgets includes
#include <qMRMLUtils.h>

// MRML includes
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLAnnotationHierarchyNode.h"
#include "vtkMRMLAnnotationFiducialNode.h"
#include "vtkMRMLFiducial.h"
#include "vtkMRMLFiducialListNode.h"
#include "vtkMRMLFiducialListStorageNode.h"
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLDisplayNode.h"
#include "vtkMRMLModelDisplayNode.h"

// VTK includes
#include <vtkCallbackCommand.h>
#include "vtkPolyData.h"
#include <vtkLineSource.h>
#include <vtkSphereSource.h>
#include <vtkCylinderSource.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include "vtkCollection.h"
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include "vtkSTLReader.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataNormals.h"
#include <vtkExodusIIReader.h>

static int count=0;
//-----------------------------------------------------------------------------
class qSlicerLITTPlanV2ModuleWidgetPrivate: public Ui_qSlicerLITTPlanV2Module
{
  Q_DECLARE_PUBLIC(qSlicerLITTPlanV2ModuleWidget);
protected:
  qSlicerLITTPlanV2ModuleWidget* const q_ptr;
public:
  qSlicerLITTPlanV2ModuleWidgetPrivate(qSlicerLITTPlanV2ModuleWidget& object);
  vtkSlicerTransformLogic*      logic()const;
  QButtonGroup*                 CoordinateReferenceButtonGroup;
  vtkMRMLLinearTransformNode*   MRMLTransformNode;
};
//-----------------------------------------------------------------------------
qSlicerLITTPlanV2ModuleWidgetPrivate::qSlicerLITTPlanV2ModuleWidgetPrivate(qSlicerLITTPlanV2ModuleWidget& object)
  : q_ptr(&object)
{
  this->CoordinateReferenceButtonGroup = 0;
  this->MRMLTransformNode = 0;
}
//-----------------------------------------------------------------------------
vtkSlicerTransformLogic* qSlicerLITTPlanV2ModuleWidgetPrivate::logic()const
{
  Q_Q(const qSlicerLITTPlanV2ModuleWidget);
  return vtkSlicerTransformLogic::SafeDownCast(q->logic());
}
//-----------------------------------------------------------------------------
qSlicerLITTPlanV2ModuleWidget::qSlicerLITTPlanV2ModuleWidget(QWidget* _parentWidget)
  : Superclass(_parentWidget)
  , d_ptr(new qSlicerLITTPlanV2ModuleWidgetPrivate(*this))
{
	SetParametersToDefault();	
	tissueIndex=0;

    this->CallBack1 = vtkSmartPointer<vtkCallbackCommand>::New();
    this->CallBack2 = vtkSmartPointer<vtkCallbackCommand>::New();

    this->fnode0 = vtkSmartPointer<vtkMRMLAnnotationFiducialNode>::New();
    this->fnode1 = vtkSmartPointer<vtkMRMLAnnotationFiducialNode>::New();
    femWatcher=new QFileSystemWatcher();

    observersAreActive=false;
    fiducialsAreValid=false;
    applicatorVisible=true;
    pathVisible=true;
    startPoint[0]=0.0; startPoint[1]=0.0; startPoint[2]=0.0;
    targetPoint[0]=0.0; targetPoint[1]=0.0; targetPoint[2]=0.0;
    prevStartPoint[0]=0.0; prevStartPoint[1]=0.0; prevStartPoint[2]=0.0;
    prevTargetPoint[0]=0.0; prevTargetPoint[1]=0.0; prevTargetPoint[2]=0.0;
}
//-----------------------------------------------------------------------------
qSlicerLITTPlanV2ModuleWidget::~qSlicerLITTPlanV2ModuleWidget()
{	
	if (this->mrmlScene())
    {
		this->mrmlScene()->RemoveAllObservers();
		//this->mrmlScene()->RemoveObserver(this->CallBack1);
		//this->mrmlScene()->RemoveObserver(this->CallBack2);
    }
	this->observersAreActive=false;	
    delete femWatcher;
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::setup()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  d->setupUi(this);

  CallBack1->SetClientData(this);
  CallBack2->SetClientData(this);
  
  // Initialize the parameters textboxes
  d->ThermalConductivity->setText("0.00");
  d->TissuePerfusion->setText("0.00");
  d->OpticalAbsorption->setText("0.00");
  d->OpticalScattering->setText("0.00");
  d->OpticalAnisotrophy->setText("0.00");
  d_ptr->ThermalDose->setText(QString::number(thermalDose));
  d_ptr->powerCoefficient->setText(QString::number(powerCoefficient));
  d_ptr->Power->setText(QString::number(power));
  d_ptr->nelemtip->setText(QString::number(nelemtip));  
  d_ptr->radius->setText(QString::number(radius));
  d_ptr->wmtLabel->setText("4");
  d_ptr->gmtLabel->setText("2");
  d_ptr->csfLabel->setText("3");
  d_ptr->tumLabel->setText("5");
  d_ptr->vesLabel->setText("1");
  d_ptr->edeLabel->setText("6");
  d_ptr->othLabel->setText("7");
  d_ptr->segFolder->setText("./Patients/Patient1.vtk");
  d_ptr->iniFolder->setText("./PatientData.ini");
  d_ptr->femFolder->setText("./fem.stl");
  d_ptr->finFolder->setText("./fem.finish");
  d_ptr->vtkFolder->setText("./fem.vtk");

  // When MRML scene changed load its contens to inputFiducials node selector (the list of the fiducials)
  connect(this, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)), d->inputFiducialsNodeSelector, SLOT(setMRMLScene(vtkMRMLScene*)));

  connect(femWatcher, SIGNAL(fileChanged(QString)), this, SLOT(LoadMesh()));

  // Load label map button is clicked
  connect(d->btnLabelmap, SIGNAL(clicked()), this, SLOT(LoadLabelMap()));

  // When the fiducials are selected load them into the class
  connect(d->inputFiducialsNodeSelector, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(onFiducalNodeChanged(vtkMRMLNode*)));
 
  // Hide Path button is clicked
  connect(d->btnHide, SIGNAL(clicked()), this, SLOT(HidePath()));

  // Show Path button is clicked
  connect(d->btnPath, SIGNAL(clicked()), this, SLOT(CreatePath()));
  
  // Load Applicator button is clicked
  connect(d->btnApplicator, SIGNAL(clicked()), this, SLOT(LoadApplicator()));
  
  //connect(d->btnFEM, SIGNAL(clicked()), this, SLOT(GetFEM()));
  connect(d->btnBurn, SIGNAL(clicked()), this, SLOT(onBtnBurnClicked()));
  connect(d->TissueSelector, SIGNAL(currentIndexChanged(int)), this, SLOT(onTissueTypeChanged()));
  connect(d->btnDefault, SIGNAL(clicked()), this, SLOT(onBtnDefaultClicked()));
  
  // Add coordinate reference button to a button group
  d->CoordinateReferenceButtonGroup = new QButtonGroup(d->CoordinateReferenceGroupBox);
  d->CoordinateReferenceButtonGroup->addButton(d->GlobalRadioButton, qMRMLTransformSliders::GLOBAL);
  d->CoordinateReferenceButtonGroup->addButton(d->LocalRadioButton, qMRMLTransformSliders::LOCAL);
  
  // Connect button group
  this->connect(d->CoordinateReferenceButtonGroup, SIGNAL(buttonPressed(int)), SLOT(onCoordinateReferenceButtonPressed(int)));

  // Connect identity button
  this->connect(d->IdentityPushButton, SIGNAL(clicked()), SLOT(identity()));

  // Connect revert button
  this->connect(d->InvertPushButton, SIGNAL(clicked()), SLOT(invert()));

  // Connect node selector with module itself
  this->connect(d->TransformNodeSelector, SIGNAL(currentNodeChanged(vtkMRMLNode*)), SLOT(onNodeSelected(vtkMRMLNode*)));

  // Connect minimum and maximum from the translation sliders to the matrix
  this->connect(d->TranslationSliders, SIGNAL(rangeChanged(double,double)), SLOT(onTranslationRangeChanged(double,double)));

  // Notify the matrix of the current translation min/max values
  this->onTranslationRangeChanged(d->TranslationSliders->minimum(), d->TranslationSliders->maximum());

  // Transform nodes connection
  this->connect(d->TransformToolButton, SIGNAL(clicked()), SLOT(transformSelectedNodes())); 
  this->connect(d->UntransformToolButton, SIGNAL(clicked()), SLOT(untransformSelectedNodes()));

  // Icons
  QIcon rightIcon = QApplication::style()->standardIcon(QStyle::SP_ArrowRight);
  d->TransformToolButton->setIcon(rightIcon);

  QIcon leftIcon = QApplication::style()->standardIcon(QStyle::SP_ArrowLeft);
  d->UntransformToolButton->setIcon(leftIcon);

  this->onNodeSelected(0);
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onFiducalNodeChanged(vtkMRMLNode* node)
{
	std::vector<vtkMRMLNode*> nodes;
	
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
				
	if (nodes.size() > 1 && (d_ptr->inputFiducialsNodeSelector->currentNode() != 0))
	{		
		for (unsigned int i=0; i<nodes.size(); i++)
		{			
			vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
			vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
			
			if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
				&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL)
			{
				fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
				fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));
				
				if(!observersAreActive)
				{
					CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
					CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);
					
					fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
					fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);
					observersAreActive=true;
				}

				fnode0->GetFiducialCoordinates(startPoint); 
				fnode1->GetFiducialCoordinates(targetPoint);				
			}
		}
	}
	else
	{
		std::vector<vtkMRMLNode*> modelNodes;
		this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:Applicator update its polydata
			{			
				//this->mrmlScene()->RemoveNode(this->mrmlScene()->GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
				//this->mrmlScene()->RemoveNode(modelNodes[j]);
				break;
			}
		}	
		this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);

		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "Applicator")==0) // There is already a model:Applicator update its polydata
			{			
				//this->mrmlScene()->RemoveNode(this->mrmlScene()->GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
				//his->mrmlScene()->RemoveNode(modelNodes[j]);
				break;
			}
		}
		modelNodes.clear();
	}
	nodes.clear();		
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadLabelMap()
{
    // Read the labelmap segmentation
    qSlicerIO::IOProperties parameters;
    parameters["fileName"] = QString(d_ptr->segFolder->text());
    parameters["labelmap"] = true;
    //parameters["center"] = true;
    parameters["discardOrientation"] = true;

    qSlicerCoreApplication::application()->coreIOManager()->loadNodes(qSlicerIO::VolumeFile, parameters);

}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::CreatePath()
{
	std::vector<vtkMRMLNode*> nodes;	
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
				
	if (nodes.size() < 2 || d_ptr->inputFiducialsNodeSelector->currentNode() == 0) // At least 2 nodes must exist(fiducials to create path)
	{
		nodes.clear(); 
		return;
	}

	for (unsigned int i=0; i<nodes.size(); i++)
	{			
		vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
		vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
		
		if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
			&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL) // There are at least two fiducial child nodes
		{
			fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
			fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));

			if(fnode0->GetFiducialCoordinates(startPoint) && fnode1->GetFiducialCoordinates(targetPoint))
			{
				CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
				CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);

				fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
				fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);

				std::vector<vtkMRMLNode*> modelNodes;
				this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
				
				// Define a line for the path
				vtkSmartPointer<vtkLineSource> line=vtkSmartPointer<vtkLineSource>::New();
				vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();						
				line->SetPoint1(startPoint[0], startPoint[1], startPoint[2]);
				line->SetPoint2(targetPoint[0], targetPoint[1], targetPoint[2]);									
				polydata=line->GetOutput();

				// Define a cylinder for the applicator
				//vtkSmartPointer<vtkCylinderSource> cylinder=vtkSmartPointer<vtkCylinderSource>::New(); // Create the applicator as a cylinder				
				//double x=targetPoint[0]-startPoint[0]; // Set the direction vector of the applicator
				//double y=targetPoint[1]-startPoint[1];
				//double z=targetPoint[2]-startPoint[2];
				//double h=50; // Height of cylinder
				//double r=2;  // Radius of cylinder
				//double angle=(180/3.1415926)*acos(y/sqrt(x*x+y*y+z*z)); // The angle between the applicator and the y-axis (default axis of creation)
				//cylinder->SetHeight(h); 
				//cylinder->SetRadius(r); 
				//cylinder->SetCenter(0.0, h/2., 0.0); 
				//cylinder->SetResolution(16);								
				
				// Define a transform for the cylinder
				vtkSmartPointer<vtkTransform> transform=vtkSmartPointer<vtkTransform>::New();
				//transform->Translate(startPoint[0], startPoint[1], startPoint[2]);
				//transform->RotateWXYZ(angle, z, 0, -1*x);				
				//vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();								
				//polydata=cylinder->GetOutput();

				vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode=vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
				transformNode->SetAndObserveMatrixTransformToParent( transform->GetMatrix() );
				transformNode->SetName("Transform Path");


				bool found=false;
									
				for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
				{
					if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:Applicator update its polydata
					{	
						this->mrmlScene()->RemoveNode(this->mrmlScene()->
							GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
						this->mrmlScene()->AddNode(transformNode);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObserveTransformNodeID(transformNode->GetID());
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);
											
						pathVisible=true;
						found=true;
						break;
					}
				}				

				if(!found)
				{
					this->mrmlScene()->AddNode(transformNode);
					vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();							
					modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
					modelNode->SetScene(this->mrmlScene());
					modelNode->SetName("ApplicatorPath");
					modelNode->SetAndObservePolyData(polydata);					

					vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();							
                    modelDisplayNode->SetColor(0,1,0); // Green
					modelDisplayNode->SetVisibility(1);
					pathVisible=true;

					this->mrmlScene()->AddNode(modelDisplayNode);
					modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

					this->mrmlScene()->AddNode(modelNode);	
				}
				modelNodes.clear();
			}
		}
	}	
	nodes.clear();	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::HidePath()
{
	std::vector<vtkMRMLNode*> modelNodes;
	this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
	for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
	{
		if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:Applicator update its polydata
		{			
			this->mrmlScene()->RemoveNode(this->mrmlScene()->GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
			this->mrmlScene()->RemoveNode(modelNodes[j]);
			break;
		}
	}	
	modelNodes.clear();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadApplicator()
{	
	/*d_ptr->btnPath->setText(QString::number(count));
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\laserApplicator.stl");*/
	
    //vtkSmartPointer<vtkExodusIIReader> reader = vtkSmartPointer<vtkExodusIIReader>::New();
	//   reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\fem.e");
	//reader->ExodusModelMetadataOn();
	//reader->UpdateInformation();
	////reader->SetPointResultArrayStatus("u", 0);
	//vtkSmartPointer<vtkContourFilter> filter = vtkSmartPointer<vtkContourFilter>::New();
	//filter->SetInputConnection(reader->GetOutputPort());
	////filter->SetValue(0, 50);
	//vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	//normals->SetInputConnection(filter->GetOutputPort());
	////normals->SetFeatureAngle(60.0);

	//vtkMRMLLinearTransformNode* transformNode=vtkMRMLLinearTransformNode::New();
	//this->mrmlScene()->AddNode(transformNode);

	//vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
 //   mat->Identity();
	//double x=targetPoint[0]-startPoint[0];
	//double y=targetPoint[1]-startPoint[1];
	//double z=targetPoint[2]-startPoint[2];
	//double h=sqrt(x*x+y*y+z*z);
	//double angle=acos(y/h);
	//angle=angle*(180/3.1415926); // Convert from radians to degrees

	////double h1=sqrt(x*x+y*y);
	////double h2=sqrt(z*z+y*y);
	////double cosT=y/h2, sinT=z/h2, cosA=y/h1, sinA=x/h1;

	///*mat->SetElement(0, 0, cosA);      mat->SetElement(0, 1, -1*sinA);   mat->SetElement(0, 2, 0);
	//mat->SetElement(1, 0, cosT*sinA); mat->SetElement(1, 1, cosT*cosA); mat->SetElement(1, 2, -1*sinT);
	//mat->SetElement(2, 0, sinT*sinA); mat->SetElement(2, 1, sinT*cosA); mat->SetElement(2, 2, cosT);		
	//mat->SetElement(0,3,startPoint[0]);
	//mat->SetElement(1,3,startPoint[1]);
	//mat->SetElement(2,3,startPoint[2]);
	//transformNode->SetAndObserveMatrixTransformToParent( mat.GetPointer() );
	//*/

	//vtkSmartPointer<vtkTransform> transform=vtkSmartPointer<vtkTransform>::New();
	//
	//transform->Translate(startPoint[0] + x/2, startPoint[1] + y/2, startPoint[2] + z/2);
	//transform->RotateWXYZ(angle, z, 0, -1*x);

	////transform->RotateX(atan(z/y)*(180/3.14));
	//d_ptr->btnBurn->setText(QString::number(angle));
	////transformNode->ApplyTransform(transform);

	////transformNode->ApplyTransformMatrix(transform->GetMatrix());
	//transformNode->SetAndObserveMatrixTransformToParent( transform->GetMatrix() );
	//
	//vtkMRMLModelNode* modelNode=vtkMRMLModelNode::New();
	//modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
	//modelNode->SetScene(this->mrmlScene());
	//modelNode->SetName("Applicator");
	//modelNode->SetAndObservePolyData(reader->GetOutput());	
 //   reader->Update();

	//vtkMRMLModelDisplayNode* modelDisplayNode=vtkMRMLModelDisplayNode::New();
	//modelDisplayNode->SetColor(1,0,0); // green
	//modelDisplayNode->SetScene(this->mrmlScene());
	//this->mrmlScene()->AddNode(modelDisplayNode);
	//modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

	//modelDisplayNode->SetPolyData(modelNode->GetPolyData());
	//this->mrmlScene()->AddNode(modelNode);

	//transformNode->Delete();
	//modelNode->Delete();
	//modelDisplayNode->Delete();	    

	std::vector<vtkMRMLNode*> nodes;	
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
				
	if (nodes.size() < 2 || d_ptr->inputFiducialsNodeSelector->currentNode() == 0) // At least 2 nodes must exist(fiducials to create path)
	{
		nodes.clear(); 
		return;
	}

	for (unsigned int i=0; i<nodes.size(); i++)
	{			
		vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
		vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
		
		if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
			&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL) // There are at least two fiducial child nodes
		{
			fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
			fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));

			if(fnode0->GetFiducialCoordinates(startPoint) && fnode1->GetFiducialCoordinates(targetPoint))
			{
				CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
				CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);

				fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
				fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);

				std::vector<vtkMRMLNode*> modelNodes;
				this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
				
				// Define a cylinder for the applicator
                vtkSmartPointer<vtkCylinderSource> cylinder=vtkSmartPointer<vtkCylinderSource>::New(); // Create the applicator tip as a cylinder
                vtkSmartPointer<vtkCylinderSource> cylinder2=vtkSmartPointer<vtkCylinderSource>::New(); // Create the applicator as a cylinder
				double x=targetPoint[0]-startPoint[0]; // Set the direction vector of the applicator
				double y=targetPoint[1]-startPoint[1];
				double z=targetPoint[2]-startPoint[2];
                double h=10, h2=sqrt(x*x+y*y+z*z); // Height of tip and applicator
                double r=1.1, r2=1;  // Radius of tip and applicator
				double angle=(180/3.1415926)*acos(y/sqrt(x*x+y*y+z*z)); // The angle between the applicator and the y-axis (default axis of creation)

                cylinder->SetHeight(h); // Define applicator tip
                cylinder->SetRadius(r);
                cylinder->SetCenter(0.0, 0.0, 0.0);
                cylinder->SetResolution(16);

                cylinder2->SetHeight(h2); // Define applicator
                cylinder2->SetRadius(r2);
                cylinder2->SetCenter(0.0, 0.0, 0.0);
                cylinder2->SetResolution(16);
				
                // Define a transform for the tip
				vtkSmartPointer<vtkTransform> transform=vtkSmartPointer<vtkTransform>::New();
                transform->Translate(startPoint[0], startPoint[1], startPoint[2]);
				transform->RotateWXYZ(angle, z, 0, -1*x);				
				vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();								
				polydata=cylinder->GetOutput();

				vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode=vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
				transformNode->SetAndObserveMatrixTransformToParent( transform->GetMatrix() );
                transformNode->SetName("Transform Applicator-tip");

                // Define a transform for the applicator
                vtkSmartPointer<vtkTransform> transform2=vtkSmartPointer<vtkTransform>::New();
                transform2->Translate((startPoint[0]+targetPoint[0])/2., (startPoint[1]+targetPoint[1])/2., (startPoint[2]+targetPoint[2])/2.);
                transform2->RotateWXYZ(angle, z, 0, -1*x);
                vtkSmartPointer<vtkPolyData> polydata2=vtkSmartPointer<vtkPolyData>::New();
                polydata2=cylinder2->GetOutput();

                vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode2=vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
                transformNode2->SetAndObserveMatrixTransformToParent( transform2->GetMatrix() );
                transformNode2->SetName("Transform Applicator-body");


                bool found1=false, found2=false;
									
				for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
				{
                    if(strcmp(modelNodes[j]->GetName(), "Applicator-tip")==0) // There is already a model:Tip update its polydata
					{	
						this->mrmlScene()->RemoveNode(this->mrmlScene()->
							GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
						this->mrmlScene()->AddNode(transformNode);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObserveTransformNodeID(transformNode->GetID());
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);
											
						applicatorVisible=true;
                        found1=true;
					}
                    if(strcmp(modelNodes[j]->GetName(), "Applicator-body")==0) // There is already a model:Applicator update its polydata
                    {
                        this->mrmlScene()->RemoveNode(this->mrmlScene()->
                            GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
                        this->mrmlScene()->AddNode(transformNode2);
                        vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata2);
                        vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObserveTransformNodeID(transformNode2->GetID());
                        vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);

                        applicatorVisible=true;
                        found2=true;
                    }
                    if(found1 && found2) break;
				}				

                if(!found1 || !found2)
				{
                    // Applicator tip
					this->mrmlScene()->AddNode(transformNode);
					vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();							
					modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
					modelNode->SetScene(this->mrmlScene());
                    modelNode->SetName("Applicator-tip");
					modelNode->SetAndObservePolyData(polydata);					

					vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();							
                    modelDisplayNode->SetColor(0,0,1); // Blue
					modelDisplayNode->SetVisibility(1);
					applicatorVisible=true;

					this->mrmlScene()->AddNode(modelDisplayNode);
					modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

					this->mrmlScene()->AddNode(modelNode);	

                    // applicator body
                    this->mrmlScene()->AddNode(transformNode2);
                    vtkSmartPointer<vtkMRMLModelNode> modelNode2=vtkSmartPointer<vtkMRMLModelNode>::New();
                    modelNode2->SetAndObserveTransformNodeID(transformNode2->GetID());
                    modelNode2->SetScene(this->mrmlScene());
                    modelNode2->SetName("Applicator-body");
                    modelNode2->SetAndObservePolyData(polydata2);

                    vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode2=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
                    modelDisplayNode2->SetColor(0,0,1); // Blue
                    modelDisplayNode2->SetVisibility(1);
                    applicatorVisible=true;

                    this->mrmlScene()->AddNode(modelDisplayNode2);
                    modelNode2->SetAndObserveDisplayNodeID(modelDisplayNode2->GetID());

                    this->mrmlScene()->AddNode(modelNode2);
				}
				modelNodes.clear();
			}
		}
	}	
	nodes.clear();	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::CreateSphereAtTarget()
{	
	std::vector<vtkMRMLNode*> nodes;	
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
				
	if (nodes.size() < 2 || d_ptr->inputFiducialsNodeSelector->currentNode() == 0) // At least 2 nodes must exist(fiducials to create path)
	{
		nodes.clear(); 
		return;
	}

	for (unsigned int i=0; i<nodes.size(); i++)
	{			
		vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
		vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
		
		if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
			&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL) // There are at least two fiducial child nodes
		{
			fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
			fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));

			if(fnode0->GetFiducialCoordinates(startPoint) && fnode1->GetFiducialCoordinates(targetPoint))
			{
				CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
				CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);

				fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
				fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);

				std::vector<vtkMRMLNode*> modelNodes;
				this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
				
				// Define a cylinder for the applicator
				vtkSmartPointer<vtkSphereSource> sphere=vtkSmartPointer<vtkSphereSource>::New(); 	
				vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();
				sphere->SetCenter(targetPoint[0], targetPoint[1], targetPoint[2]);
				sphere->SetRadius(5); 	
				
				sphere->SetThetaResolution(32);
				sphere->SetPhiResolution(32);
				sphere->SetStartTheta(0.0);
				sphere->SetEndTheta(360.0);
				sphere->SetStartPhi(0.0);
				sphere->SetEndPhi(360.0);
				sphere->SetLatLongTessellation(1);
				polydata=sphere->GetOutput();

                //sleep(500);
				
				// Define a transform for the cylinder
				vtkSmartPointer<vtkTransform> transform=vtkSmartPointer<vtkTransform>::New();
				
				vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode=vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
				transformNode->SetAndObserveMatrixTransformToParent( transform->GetMatrix() );
				transformNode->SetName("Transform Treatment");


				bool found=false;
									
				for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
				{
					if(strcmp(modelNodes[j]->GetName(), "Treatment")==0) // There is already a model:Applicator update its polydata
					{	
						this->mrmlScene()->RemoveNode(this->mrmlScene()->
							GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
						this->mrmlScene()->AddNode(transformNode);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObserveTransformNodeID(transformNode->GetID());
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);
											
						//pathVisible=true;
						found=true;
						break;
					}
				}				

				if(!found)
				{
					this->mrmlScene()->AddNode(transformNode);
					vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();							
					modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
					modelNode->SetScene(this->mrmlScene());
					modelNode->SetName("Treatment");
					modelNode->SetAndObservePolyData(polydata);					

					vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();							
					modelDisplayNode->SetColor(0,0,1); // Red
					modelDisplayNode->SetVisibility(1);
					//pathVisible=true;

					this->mrmlScene()->AddNode(modelDisplayNode);
					modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

					this->mrmlScene()->AddNode(modelNode);	
				}
				modelNodes.clear();
			}
		}
	}	
	nodes.clear();	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved(vtkObject* vtk_obj, unsigned long event, void* client_data, void* call_data)
{
	vtkMRMLAnnotationFiducialNode* fnode = reinterpret_cast<vtkMRMLAnnotationFiducialNode*>(vtk_obj);
	qSlicerLITTPlanV2ModuleWidget* thisClass = reinterpret_cast<qSlicerLITTPlanV2ModuleWidget*>(client_data);
	
	//if(thisClass->applicatorVisible || thisClass->pathVisible)
	{
		fnode->GetFiducialCoordinates(thisClass->startPoint);

		/*vtkSmartPointer<vtkLineSource> line=vtkSmartPointer<vtkLineSource>::New();
		vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();						
		line->SetPoint1(thisClass->startPoint[0], thisClass->startPoint[1], thisClass->startPoint[2]);
		line->SetPoint2(thisClass->targetPoint[0], thisClass->targetPoint[1], thisClass->targetPoint[2]);									
		polydata=line->GetOutput();*/

		std::vector<vtkMRMLNode*> modelNodes;	
		thisClass->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
		
		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->applicatorVisible=false;
				break;
			}
		}

		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
            if(strcmp(modelNodes[j]->GetName(), "Applicator-tip")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
                vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetColor(0,0,1);
				//thisClass->pathVisible=false;
				break;
			}
		}
        for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
        {
            if(strcmp(modelNodes[j]->GetName(), "Applicator-body")==0) // There is already a model:ApplicatorPath update its polydata
            {
                //if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
                //vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
                //thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
                vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
                //thisClass->pathVisible=false;
                break;
            }
        }
		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "Treatment")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->pathVisible=false;
				break;
			}
		}
		modelNodes.clear();
	}

	/*for(unsigned int c=0; c<3; c++)
	{
		thisClass->prevStartPoint[c]=thisClass->startPoint[c];
	}*/

	//thisClass->CreatePath();
	
	//vtkMRMLAnnotationFiducialNode* fnode0 = reinterpret_cast<vtkMRMLAnnotationFiducialNode*>(vtk_obj);

	//if(fnode0->GetFiducialCoordinates(qSlicerLITTPlanV2ModuleWidget::startPoint))
	//{
	//}

	//vtkMRMLAnnotationFiducialNode *fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(vtk_obj);
	//vtkMRMLAnnotationFiducialNode *fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(vtk_obj);
	
	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved(vtkObject* vtk_obj, unsigned long event, void* client_data, void* call_data)
{
	vtkMRMLAnnotationFiducialNode* fnode = reinterpret_cast<vtkMRMLAnnotationFiducialNode*>(vtk_obj);
	qSlicerLITTPlanV2ModuleWidget* thisClass = reinterpret_cast<qSlicerLITTPlanV2ModuleWidget*>(client_data);
	
	//if(thisClass->applicatorVisible || thisClass->pathVisible)
	{
		fnode->GetFiducialCoordinates(thisClass->targetPoint);

		//vtkSmartPointer<vtkLineSource> line=vtkSmartPointer<vtkLineSource>::New();
		//vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();						
		//line->SetPoint1(thisClass->startPoint[0], thisClass->startPoint[1], thisClass->startPoint[2]);
		//line->SetPoint2(thisClass->targetPoint[0], thisClass->targetPoint[1], thisClass->targetPoint[2]);									
		//polydata=line->GetOutput();

		std::vector<vtkMRMLNode*> modelNodes;	
		thisClass->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
		
		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->applicatorVisible=false;
				break;
			}
		}

		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
            if(strcmp(modelNodes[j]->GetName(), "Applicator-tip")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
                vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetColor(0,0,1);
				//thisClass->pathVisible=false;
				break;
			}
		}
        for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
        {
            if(strcmp(modelNodes[j]->GetName(), "Applicator-body")==0) // There is already a model:ApplicatorPath update its polydata
            {
                //if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
                //vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
                //thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
                vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
                //thisClass->pathVisible=false;
                break;
            }
        }

		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "Treatment")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->pathVisible=false;
				break;
			}
		}
		modelNodes.clear();
	}

	/*for(unsigned int c=0; c<3; c++)
	{
		thisClass->prevStartPoint[c]=thisClass->startPoint[c];
	}*/

	//thisClass->CreatePath();
	
	//vtkMRMLAnnotationFiducialNode* fnode0 = reinterpret_cast<vtkMRMLAnnotationFiducialNode*>(vtk_obj);

	//if(fnode0->GetFiducialCoordinates(qSlicerLITTPlanV2ModuleWidget::startPoint))
	//{
	//}

	//vtkMRMLAnnotationFiducialNode *fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(vtk_obj);
	//vtkMRMLAnnotationFiducialNode *fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(vtk_obj);
	
	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onTissueTypeChanged()
{
    LoadParametersFromGUI(); // First get the parameter values from screen (the values before the new tissue selection)
    LoadParametersToGUI();	 // Now load the selected tissue`s parameters to the screen
    tissueIndex=d_ptr->TissueSelector->currentIndex();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::SetParametersToDefault()
{
	this->wmtParameters[0]=0.502;
	this->wmtParameters[1]=6.57;
	this->wmtParameters[2]=500;
	this->wmtParameters[3]=14000;
	this->wmtParameters[4]=0.88;

	this->gmtParameters[0]=0.527;
	this->gmtParameters[1]=6.57;
	this->gmtParameters[2]=500;
	this->gmtParameters[3]=14000;
	this->gmtParameters[4]=0.88;

	this->csfParameters[0]=0.624;
	this->csfParameters[1]=1.0;
	this->csfParameters[2]=12;
	this->csfParameters[3]=0.01;
	this->csfParameters[4]=0.9;

	this->tumParameters[0]=0.53;
	this->tumParameters[1]=2.0;
	this->tumParameters[2]=51;
	this->tumParameters[3]=8000;
	this->tumParameters[4]=0.95;	

    this->vesParameters[0]=0.7;
    this->vesParameters[1]=3.0;
    this->vesParameters[2]=45;
    this->vesParameters[3]=9000;
    this->vesParameters[4]=0.85;

    this->edeParameters[0]=0.6;
    this->edeParameters[1]=4.0;
    this->edeParameters[2]=49;
    this->edeParameters[3]=8500;
    this->edeParameters[4]=0.9;

    this->othParameters[0]=0.60;
    this->othParameters[1]=1.0;
    this->othParameters[2]=30;
    this->othParameters[3]=7000;
    this->othParameters[4]=0.50;

    power=10.0;
    thermalDose=42.0;
    powerCoefficient=0.50;
    nelemtip=10;
    radius=0.75;
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadParametersToGUI()
{
	switch(d_ptr->TissueSelector->currentIndex())
	{
    case 0: //no tissue
		d_ptr->ThermalConductivity->setText("0.00");
		d_ptr->TissuePerfusion->setText("0.00");
		d_ptr->OpticalAbsorption->setText("0.00");
		d_ptr->OpticalScattering->setText("0.00");
		d_ptr->OpticalAnisotrophy->setText("0.00"); 
		break;
    case 1: // white
		d_ptr->ThermalConductivity->setText(QString::number(wmtParameters[0]));
		d_ptr->TissuePerfusion->setText(QString::number(wmtParameters[1]));
		d_ptr->OpticalAbsorption->setText(QString::number(wmtParameters[2]));
		d_ptr->OpticalScattering->setText(QString::number(wmtParameters[3]));
		d_ptr->OpticalAnisotrophy->setText(QString::number(wmtParameters[4])); 
		break;
    case 2: // grey
		d_ptr->ThermalConductivity->setText(QString::number(gmtParameters[0]));
		d_ptr->TissuePerfusion->setText(QString::number(gmtParameters[1]));
		d_ptr->OpticalAbsorption->setText(QString::number(gmtParameters[2]));
		d_ptr->OpticalScattering->setText(QString::number(gmtParameters[3]));
		d_ptr->OpticalAnisotrophy->setText(QString::number(gmtParameters[4])); 
		break;
    case 3: // csf
		d_ptr->ThermalConductivity->setText(QString::number(csfParameters[0]));
		d_ptr->TissuePerfusion->setText(QString::number(csfParameters[1]));
		d_ptr->OpticalAbsorption->setText(QString::number(csfParameters[2]));
		d_ptr->OpticalScattering->setText(QString::number(csfParameters[3]));
		d_ptr->OpticalAnisotrophy->setText(QString::number(csfParameters[4])); 
		break;
    case 4: // tumor
		d_ptr->ThermalConductivity->setText(QString::number(tumParameters[0]));
		d_ptr->TissuePerfusion->setText(QString::number(tumParameters[1]));
		d_ptr->OpticalAbsorption->setText(QString::number(tumParameters[2]));
		d_ptr->OpticalScattering->setText(QString::number(tumParameters[3]));
		d_ptr->OpticalAnisotrophy->setText(QString::number(tumParameters[4]));  
		break;
    case 5: // vessel
        d_ptr->ThermalConductivity->setText(QString::number(vesParameters[0]));
        d_ptr->TissuePerfusion->setText(QString::number(vesParameters[1]));
        d_ptr->OpticalAbsorption->setText(QString::number(vesParameters[2]));
        d_ptr->OpticalScattering->setText(QString::number(vesParameters[3]));
        d_ptr->OpticalAnisotrophy->setText(QString::number(vesParameters[4]));
        break;
    case 6: // edema
        d_ptr->ThermalConductivity->setText(QString::number(edeParameters[0]));
        d_ptr->TissuePerfusion->setText(QString::number(edeParameters[1]));
        d_ptr->OpticalAbsorption->setText(QString::number(edeParameters[2]));
        d_ptr->OpticalScattering->setText(QString::number(edeParameters[3]));
        d_ptr->OpticalAnisotrophy->setText(QString::number(edeParameters[4]));
        break;
    case 7: // other
        d_ptr->ThermalConductivity->setText(QString::number(othParameters[0]));
        d_ptr->TissuePerfusion->setText(QString::number(othParameters[1]));
        d_ptr->OpticalAbsorption->setText(QString::number(othParameters[2]));
        d_ptr->OpticalScattering->setText(QString::number(othParameters[3]));
        d_ptr->OpticalAnisotrophy->setText(QString::number(othParameters[4]));
        break;
	default:
		break;
	}
	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadParametersFromGUI()
{
	switch(tissueIndex)
	{
	case 0:
		break;
	case 1:
		this->wmtParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
		this->wmtParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
		this->wmtParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
		this->wmtParameters[3]=d_ptr->OpticalScattering->text().toDouble();
		this->wmtParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
		break;
	case 2:
		this->gmtParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
		this->gmtParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
		this->gmtParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
		this->gmtParameters[3]=d_ptr->OpticalScattering->text().toDouble();
		this->gmtParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
		break;
	case 3:
		this->csfParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
		this->csfParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
		this->csfParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
		this->csfParameters[3]=d_ptr->OpticalScattering->text().toDouble();
		this->csfParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
		break;

    case 5:
        this->vesParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
        this->vesParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
        this->vesParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
        this->vesParameters[3]=d_ptr->OpticalScattering->text().toDouble();
        this->vesParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
        break;
    case 6:
        this->edeParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
        this->edeParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
        this->edeParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
        this->edeParameters[3]=d_ptr->OpticalScattering->text().toDouble();
        this->edeParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
        break;
    case 7:
        this->othParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
        this->othParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
        this->othParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
        this->othParameters[3]=d_ptr->OpticalScattering->text().toDouble();
        this->othParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
        break;
	default:
		break;
	}
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onBtnDefaultClicked()
{
	SetParametersToDefault();
	LoadParametersToGUI();
	d_ptr->ThermalDose->setText(QString::number(thermalDose));
	d_ptr->Power->setText(QString::number(power));
    /*d_ptr->nelemtip->setText(QString::number(nelemtip));
    d_ptr->radius->setText(QString::number(radius));
    d_ptr->wmtLabel->setText("1");
    d_ptr->gmtLabel->setText("2");
    d_ptr->csfLabel->setText("3");
    d_ptr->tumLabel->setText("4");
    d_ptr->vesLabel->setText("5");
    d_ptr->edeLabel->setText("6");
    d_ptr->othLabel->setText("7");*/
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onBtnApplicatorClicked()
{
	d_ptr->btnApplicator->setText("alo");
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\laserApplicator.stl");
	
	//vtkSmartPointer<vtkExodusIIReader> reader = vtkSmartPointer<vtkExodusIIReader>::New();
 //   reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\fem.e");
	//reader->ExodusModelMetadataOn();
	//reader->UpdateInformation();
	////reader->SetPointResultArrayStatus("u", 0);
	//vtkSmartPointer<vtkContourFilter> filter = vtkSmartPointer<vtkContourFilter>::New();
	//filter->SetInputConnection(reader->GetOutputPort());
	////filter->SetValue(0, 50);
	//vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	//normals->SetInputConnection(filter->GetOutputPort());
	////normals->SetFeatureAngle(60.0);

	vtkMRMLLinearTransformNode* transformNode=vtkMRMLLinearTransformNode::New();
	this->mrmlScene()->AddNode(transformNode);

	vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
    mat->Identity();
	//mat->SetElement(2,3,86.0);
    /*for( size_t p = 0; p < 4; p++ )
      {
      vtkFloatingPointType point = normal[p];
      mat->SetElement(static_cast<int>(p), 0, (sign * point) );
      }*/
    /*vtkFloatingPointType oneAndAlpha = 1.0 + mat->GetElement(0, 0);
    mat->SetElement(0, 1, -1.0 * mat->GetElement(1, 0) );
    mat->SetElement(0, 2, (-1.0 * (mat->GetElement(2, 0) ) ) );
    mat->SetElement(2, 1, (-1.0 * (mat->GetElement(1, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );
    mat->SetElement(1, 2, (-1.0 * (mat->GetElement(1, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );
    mat->SetElement(1, 1, (1.0  - (mat->GetElement(1, 0) * (mat->GetElement(1, 0) / oneAndAlpha) ) ) );
    mat->SetElement(2, 2, (1.0  - (mat->GetElement(2, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );*/


	transformNode->SetAndObserveMatrixTransformToParent( mat.GetPointer() );

	vtkMRMLModelNode* modelNode=vtkMRMLModelNode::New();
	modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
	modelNode->SetScene(this->mrmlScene());
	modelNode->SetName("Applicator");
	modelNode->SetAndObservePolyData(reader->GetOutput());	
    reader->Update();

	vtkMRMLModelDisplayNode* modelDisplayNode=vtkMRMLModelDisplayNode::New();
	modelDisplayNode->SetColor(1,0,0); // green
	modelDisplayNode->SetScene(this->mrmlScene());
	this->mrmlScene()->AddNode(modelDisplayNode);
	modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

	modelDisplayNode->SetPolyData(modelNode->GetPolyData());
	this->mrmlScene()->AddNode(modelNode);

	transformNode->Delete();
	modelNode->Delete();
	modelDisplayNode->Delete();	    
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onBtnBurnClicked() // Export the created points here!! FEM will be read too..
{
    // --- Get all the model nodes. Find the ``Applicator-tip`` and change its color from blue to red ---
    std::vector<vtkMRMLNode*> modelNodes;
    this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
    for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
    {
        if(strcmp(modelNodes[j]->GetName(), "Applicator-tip")==0) // There is already a model:Tip update its polydata
        {
            vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetColor(1, 0, 0);
            break;
        }
    }
    modelNodes.clear();

    // --- Update all the parameters from GUI ---
	LoadParametersFromGUI(); // Tissue parameters
	this->thermalDose=d_ptr->ThermalDose->text().toDouble(); // Global parameter: Thermal dose	
    this->nelemtip=d_ptr->nelemtip->text().toDouble(); // Global parameter: number of elements
    this->powerCoefficient=d_ptr->powerCoefficient->text().toDouble(); // Global parameter: power coefficient
    this->power=(d_ptr->Power->text().toDouble())*this->powerCoefficient; // Global parameter: power
    this->radius=d_ptr->radius->text().toDouble(); // Global parameter: power
		
	// Write all the parameter values and the two guiding points to .ini file
    ofstream outputFile; // Define a text file

    QByteArray byteArray = d_ptr->iniFolder->text().toUtf8();
    const char* cString = byteArray.constData();

    char folder[400]="";
    strcat(folder, cString); // Fill the folder with cString which is d_ptr->iniFolder->text()
    outputFile.open(folder); // Create the text file

	
    outputFile<<"# needs to be in meters"<<endl;
    outputFile<<"# (x_0,y_0,z_0) is applicator centroid"<<endl;
    outputFile<<"# (x_1,y_1,z_1) is applicator tip"<<endl;

    outputFile<<"[probe]"<<endl;
    outputFile<<"x_0 = "<<startPoint[0]/1000.<<endl;
    outputFile<<"y_0 = "<<startPoint[1]/1000.<<endl;
    outputFile<<"z_0 = "<<startPoint[2]/1000.<<endl;
    outputFile<<"x_1 = "<<targetPoint[0]/1000.<<endl;
    outputFile<<"y_1 = "<<targetPoint[1]/1000.<<endl;
    outputFile<<"z_1 = "<<targetPoint[2]/1000.<<endl;

    outputFile<<"nelemtip = "<<nelemtip<<endl;
    outputFile<<"radius = "<<radius/1000.<<endl;

    outputFile<<"[timestep]"<<endl;
    outputFile<<"power="<<power<<endl;
    outputFile<<"thermalDose="<<thermalDose<<endl;

    outputFile<<"[thermal_conductivity]"<<endl;
    outputFile<<"k_0_healthy="<<othParameters[0]<<endl;
    outputFile<<"k_0_grey="<<gmtParameters[0]<<endl;
    outputFile<<"k_0_white="<<wmtParameters[0]<<endl;
    outputFile<<"k_0_csf="<<csfParameters[0]<<endl;
    outputFile<<"k_0_tumor="<<tumParameters[0]<<endl;
    outputFile<<"k_0_vessel="<<vesParameters[0]<<endl;
    outputFile<<"k_0_edema="<<edeParameters[0]<<endl;

    outputFile<<"[perfusion]"<<endl;
    outputFile<<"w_0_healthy="<<othParameters[1]<<endl;
    outputFile<<"w_0_grey="<<gmtParameters[1]<<endl;
    outputFile<<"w_0_white="<<wmtParameters[1]<<endl;
    outputFile<<"w_0_csf="<<csfParameters[1]<<endl;
    outputFile<<"w_0_tumor="<<tumParameters[1]<<endl;
    outputFile<<"w_0_vessel="<<vesParameters[1]<<endl;
    outputFile<<"w_0_edema="<<edeParameters[1]<<endl;

    outputFile<<"[optical]"<<endl;
    outputFile<<"mu_a_healthy="<<othParameters[2]<<endl;
    outputFile<<"mu_a_grey="<<gmtParameters[2]<<endl;
    outputFile<<"mu_a_white="<<wmtParameters[2]<<endl;
    outputFile<<"mu_a_csf="<<csfParameters[2]<<endl;
    outputFile<<"mu_a_tumor="<<tumParameters[2]<<endl;
    outputFile<<"mu_a_vessel="<<vesParameters[2]<<endl;
    outputFile<<"mu_a_edema="<<edeParameters[2]<<endl;

    outputFile<<"mu_s_healthy="<<othParameters[3]<<endl;
    outputFile<<"mu_s_grey="<<gmtParameters[3]<<endl;
    outputFile<<"mu_s_white="<<wmtParameters[3]<<endl;
    outputFile<<"mu_s_csf="<<csfParameters[3]<<endl;
    outputFile<<"mu_s_tumor="<<tumParameters[3]<<endl;
    outputFile<<"mu_s_vessel="<<vesParameters[3]<<endl;
    outputFile<<"mu_s_edema="<<edeParameters[3]<<endl;

    outputFile<<"anfact_healthy="<<othParameters[4]<<endl;
    outputFile<<"anfact_grey="<<gmtParameters[4]<<endl;
    outputFile<<"anfact_white="<<wmtParameters[4]<<endl;
    outputFile<<"anfact_csf="<<csfParameters[4]<<endl;
    outputFile<<"anfact_tumor="<<tumParameters[4]<<endl;
    outputFile<<"anfact_vessel="<<vesParameters[4]<<endl;
    outputFile<<"anfact_edema="<<edeParameters[4]<<endl;

    QByteArray byteArray1 = d_ptr->vesLabel->text().toUtf8();
    const char* cString1 = byteArray1.constData();

    QByteArray byteArray2 = d_ptr->wmtLabel->text().toUtf8();
    const char* cString2 = byteArray2.constData();

    QByteArray byteArray3 = d_ptr->gmtLabel->text().toUtf8();
    const char* cString3 = byteArray3.constData();

    QByteArray byteArray4 = d_ptr->csfLabel->text().toUtf8();
    const char* cString4 = byteArray4.constData();

    QByteArray byteArray5 = d_ptr->tumLabel->text().toUtf8();
    const char* cString5 = byteArray5.constData();

    QByteArray byteArray6 = d_ptr->othLabel->text().toUtf8();
    const char* cString6 = byteArray6.constData();

    QByteArray byteArray7 = d_ptr->edeLabel->text().toUtf8();
    const char* cString7 = byteArray7.constData();

    QByteArray byteArray8 = d_ptr->segFolder->text().toUtf8();
    const char* cString8 = byteArray8.constData();

    outputFile<<"[labels]"<<endl;


    outputFile<<"vessel="<<cString1<<endl;
    outputFile<<"whitematter="<<cString2<<endl;
    outputFile<<"greymatter="<<cString3<<endl;
    outputFile<<"csf="<<cString4<<endl;
    outputFile<<"tumor="<<cString5<<endl;
    outputFile<<"edema="<<cString7<<endl;
    outputFile<<"other="<<cString6<<endl;

    outputFile<<"[initial_condition]"<<endl;
    outputFile<<"u_init"<<"=37.0"<<endl;


    outputFile<<"[exec]"<<endl;
    outputFile<<"segment_file ="<<cString8<<endl;
    outputFile<<"roi = [(60,250),(40,300),(10,90)]"<<endl;
    outputFile<<"contours = ["<<thermalDose<<"]"<<endl;
    outputFile<<"subsample = [3,3,2]"<<endl;
    outputFile<<";subsample = [5,5,5]"<<endl;

    outputFile.close();

    femWatcher->blockSignals((false));
    if(femWatcher->files().count()>0)  femWatcher->removePaths(femWatcher->files());
    femWatcher->addPath(d_ptr->finFolder->text());


    //LoadMesh();
    //CreateSphereAtTarget();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadMesh()
{
    count++;
    //d_ptr->ThermalDose->setText(QString::number(count));
    // Read the mesh file
    QByteArray byteArray = d_ptr->femFolder->text().toUtf8();
    const char* cString = byteArray.constData();
    char folder[400]="";
    strcat(folder, cString);

    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(folder);

    vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();
    modelNode->SetScene(this->mrmlScene());

    QByteArray byteArrayTD = d_ptr->ThermalDose->text().toUtf8();
    const char* cStringTD = byteArrayTD.constData();
    char TD[10]="";
    strcpy(TD,"Treatment_");
    strcat(TD, cStringTD);

    modelNode->SetName(TD);
    modelNode->SetAndObservePolyData(reader->GetOutput());
    reader->Update();

    vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
    if(EqualD(thermalDose, 0.0)) thermalDose=1; // no zero division below color..
    if(thermalDose>90) modelDisplayNode->SetColor(1,0,0);
    else if(thermalDose<35) modelDisplayNode->SetColor(0,0,1);
    else modelDisplayNode->SetColor(thermalDose/90.,0,35./thermalDose);
    this->mrmlScene()->AddNode(modelDisplayNode);
    modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

    this->mrmlScene()->AddNode(modelNode);

    // Read the labelmap vtk file
    //qSlicerIO::IOProperties parameters;
    //parameters["fileName"] = QString(d_ptr->vtkFolder->text());
    //parameters["labelmap"] = true;

    //qSlicerCoreApplication::application()->coreIOManager()->loadNodes(qSlicerIO::VolumeFile, parameters);

    //d_ptr->btnBurn->setText();

    femWatcher->blockSignals((true));
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onCoordinateReferenceButtonPressed(int id)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  
  qMRMLTransformSliders::CoordinateReferenceType ref =
    (id == qMRMLTransformSliders::GLOBAL) ? qMRMLTransformSliders::GLOBAL : qMRMLTransformSliders::LOCAL;
  d->TranslationSliders->setCoordinateReference(ref);
  d->RotationSliders->setCoordinateReference(ref);
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onNodeSelected(vtkMRMLNode* node)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  
  vtkMRMLLinearTransformNode* transformNode = vtkMRMLLinearTransformNode::SafeDownCast(node);

  // Enable/Disable CoordinateReference, identity buttons, MatrixViewGroupBox,
  // Min/Max translation inputs
  d->CoordinateReferenceGroupBox->setEnabled(transformNode != 0);
  d->IdentityPushButton->setEnabled(transformNode != 0);
  d->InvertPushButton->setEnabled(transformNode != 0);
  d->MatrixViewGroupBox->setEnabled(transformNode != 0);

  // Listen for Transform node changes
  this->qvtkReconnect(d->MRMLTransformNode, transformNode,
    vtkMRMLTransformableNode::TransformModifiedEvent,
    this, SLOT(onMRMLTransformNodeModified(vtkObject*)));

  QStringList nodeTypes;
  // If no transform node, it would show the entire scene, lets shown none
  // instead.
  if (transformNode == 0)
    {
    nodeTypes << QString("vtkMRMLNotANode");
    }
  d->TransformedTreeView->setNodeTypes(nodeTypes);

  // Filter the current node in the transformed tree view
  d->TransformedTreeView->setRootNode(transformNode);

  // Hide the current node in the transformable tree view
  QStringList hiddenNodeIDs;
  if (transformNode)
    {
    hiddenNodeIDs << QString(transformNode->GetID());
    }
  d->TransformableTreeView->sortFilterProxyModel()
    ->setHiddenNodeIDs(hiddenNodeIDs);
  d->MRMLTransformNode = transformNode;
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::identity()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);

  if (!d->MRMLTransformNode)
    {
    return;
    }

  d->RotationSliders->resetUnactiveSliders();
  d->MRMLTransformNode->GetMatrixTransformToParent()->Identity();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::invert()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  
  if (!d->MRMLTransformNode) { return; }

  d->RotationSliders->resetUnactiveSliders();
  d->MRMLTransformNode->GetMatrixTransformToParent()->Invert();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onMRMLTransformNodeModified(vtkObject* caller)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  
  vtkMRMLLinearTransformNode* transformNode = vtkMRMLLinearTransformNode::SafeDownCast(caller);
  if (!transformNode) { return; }

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  qMRMLUtils::getTransformInCoordinateSystem(d->MRMLTransformNode,
    this->coordinateReference() == qMRMLTransformSliders::GLOBAL, transform);

  // The matrix can be changed externally. The min/max values shall be updated 
  //accordingly to the new matrix if needed.
  vtkMatrix4x4 * mat = transform->GetMatrix();
  double min = 0.;
  double max = 0.;
  this->extractMinMaxTranslationValue(mat, min, max);
  if (min < d->TranslationSliders->minimum())
    {
    min = min - 0.3 * fabs(min);
    d->TranslationSliders->setMinimum(min);
    }
  if (max > d->TranslationSliders->maximum())
    {
    max = max + 0.3 * fabs(max);
    d->TranslationSliders->setMaximum(max);
    }
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::extractMinMaxTranslationValue(vtkMatrix4x4 * mat, double& min, double& max)
{
  if (!mat)
    {
    Q_ASSERT(mat);
    return;
    }
  for (int i=0; i <3; i++)
    {
    min = qMin(min, mat->GetElement(i,3));
    max = qMax(max, mat->GetElement(i,3));
    }
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onTranslationRangeChanged(double newMin, double newMax)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  d->MatrixWidget->setRange(newMin, newMax);
}
//-----------------------------------------------------------------------------
int qSlicerLITTPlanV2ModuleWidget::coordinateReference()const
{
  Q_D(const qSlicerLITTPlanV2ModuleWidget);
  return d->CoordinateReferenceButtonGroup->checkedId();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::setMRMLScene(vtkMRMLScene* scene)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  this->Superclass::setMRMLScene(scene);
  // If the root index is set before the scene, it will show the scene as
  // top-level item. Setting the root index to be the scene makes the nodes
  // top-level, and this can only be done once the scene is set.
  d->TransformableTreeView->setRootIndex(
    d->TransformableTreeView->sortFilterProxyModel()->mrmlSceneIndex());
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::transformSelectedNodes()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  QModelIndexList selectedIndexes =
    d->TransformableTreeView->selectionModel()->selectedRows();
  selectedIndexes = qMRMLTreeView::removeChildren(selectedIndexes);
  foreach(QModelIndex selectedIndex, selectedIndexes)
    {
    vtkMRMLTransformableNode* node = vtkMRMLTransformableNode::SafeDownCast(
    d->TransformableTreeView->sortFilterProxyModel()->
      mrmlNodeFromIndex( selectedIndex ));
    Q_ASSERT(node);
    node->SetAndObserveTransformNodeID(d->MRMLTransformNode->GetID());
    }
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::untransformSelectedNodes()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  QModelIndexList selectedIndexes =
    d->TransformedTreeView->selectionModel()->selectedRows();
  selectedIndexes = qMRMLTreeView::removeChildren(selectedIndexes);
  foreach(QModelIndex selectedIndex, selectedIndexes)
    {
    vtkMRMLTransformableNode* node = vtkMRMLTransformableNode::SafeDownCast(
    d->TransformedTreeView->sortFilterProxyModel()->
      mrmlNodeFromIndex( selectedIndex ));
    Q_ASSERT(node);
    node->SetAndObserveTransformNodeID(0);
    }
}
//-----------------------------------------------------------------------------
bool qSlicerLITTPlanV2ModuleWidget::EqualD(double a, double b)
{
	double dTolerance=0.000001;

	if(fabs(a-b)<dTolerance) 
		return true;
	else 
		return false;
}
//-----------------------------------------------------------------------------
bool qSlicerLITTPlanV2ModuleWidget::EqualP(double p1[], double p2[])
{
	if( EqualD(p1[0], p2[0]) && EqualD(p1[1], p2[1])&& EqualD(p1[2], p2[2]))
		return true;
	else 
		return false;
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::CreatePathOld()
{	
	// Alternate Approach
	//vtkMRMLNode* fiducialsNode = d_ptr->inputFiducialsNodeSelector->currentNode();

	
	//vtkMRMLFiducialListNode *fnode=vtkMRMLFiducialListNode::SafeDownCast(fiducialsNode); 
	
	//vtkMRMLAnnotationHierarchyNode *anode=vtkMRMLAnnotationHierarchyNode::SafeDownCast(fiducialsNode);

	//vtkCollection *cnodes = vtkCollection::New();	

	////if(anode->GetClassName()== "vtkMRMLAnnotationHierarchyNode")
	////{
	////	//d_ptr->label->clear();
	//	anode->GetChildrenDisplayableNodes(cnodes);

	////	
	//	int n=cnodes->GetNumberOfItems();
	//	if (n >0)
	//	{

	//		for(int i=0; i<n; i++)
	//		{
	//			d_ptr->pushButton->setText(QString::number(n));
	//			

	//			vtkMRMLFiducial *fid = vtkMRMLFiducial::SafeDownCast(cnodes->GetItemAsObject(i));
	//			if (fid != NULL)
	//			{
	//				//d_ptr->label->clear();
	//				float* xyz = fid->GetXYZ();
	//				//d_ptr->pushButton->setText(QString::number(xyz[0]));
	//			}
	//			//float coords[];
	//			//fid->GetXYZ(coords);
	//			
	//			/*vtkMRMLFiducial *fid2 = vtkMRMLFiducial::SafeDownCast(cnodes->GetItemAsObject(1));
	//			float coords2[3];
	//			fid1->GetXYZ(coords2);

	//			
	//			QString qnum= QString::number(n);
	//			d_ptr->pushButton->setText(QString::number(coords1[0]));*/
	//		}
	//	}
	////}
	//cnodes->RemoveAllItems();
	//cnodes->Delete();
	////
	// // OR Get nodes from the mrml scene...
	  
	std::vector<vtkMRMLNode*> nodes;
	std::vector<vtkMRMLNode*> modelNodes;

	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
	this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
			
	if (nodes.size() > 1 && (d_ptr->inputFiducialsNodeSelector->currentNode() != 0))
	{		
		for (unsigned int i=0; i<nodes.size(); i++)
		{			
			vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
			vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
			
			if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
				&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL)
			{
				fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
				fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));
				
				if(!observersAreActive)
				{
				
					CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
					CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);
				
					fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
					fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);
					observersAreActive=true;
				}
				
				if(fnode0->GetFiducialCoordinates(startPoint) && fnode1->GetFiducialCoordinates(targetPoint))
				{
					fiducialsAreValid=true;
					//if( (!EqualP(startPoint, prevStartPoint)) || (!EqualP(targetPoint, prevTargetPoint))) // If any of the points has been changed (=different than the previous value)
					{						
						vtkSmartPointer<vtkLineSource> line=vtkSmartPointer<vtkLineSource>::New();
						vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();						
						line->SetPoint1(startPoint[0], startPoint[1], startPoint[2]);
						line->SetPoint2(targetPoint[0], targetPoint[1], targetPoint[2]);									
						polydata=line->GetOutput();

						bool found=false;
						
						for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
						{
							if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:ApplicatorPath update its polydata
							{
								vtkMRMLModelNode* pathModel=vtkMRMLModelNode::SafeDownCast(modelNodes[j]);
								pathModel->SetAndObservePolyData(polydata);
								pathModel->GetModelDisplayNode()->SetVisibility(1);
								pathVisible=true;
								//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
								//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);
								found=true;
								break;
							}
						}
						
						if(!found) // create a new one
						{
							vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();							
							modelNode->SetScene(this->mrmlScene());
							modelNode->SetName("ApplicatorPath");							
							modelNode->SetAndObservePolyData(polydata);	
												
							vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();							
							modelDisplayNode->SetColor(1,0,0); // Red
							pathVisible=true;
							//modelDisplayNode->SetVisibility(1);

							this->mrmlScene()->AddNode(modelDisplayNode);
							modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

							this->mrmlScene()->AddNode(modelNode);													
						}				
						
						/*for(int c=0; c<3; c++)
						{
							prevStartPoint[c]=startPoint[c];
							prevTargetPoint[c]=targetPoint[c];
						}*/
					}
				}															
				break;
			}
		}
	}
	nodes.clear();	
	modelNodes.clear();
	//d_ptr->btnPath->setText(QString::number(pathVisible));
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::sleep(unsigned int mseconds)
{
    clock_t goal = mseconds + clock();
    while (goal > clock());
}
