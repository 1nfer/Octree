#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>

#include <vtkNamedColors.h>

// The ordering of the corner points on each face.
std::array<std::array<vtkIdType, 4>, 6> ordering = {{{{0, 3, 2, 1}},
                                                    {{4, 5, 6, 7}},
                                                    {{0, 1, 5, 4}},
                                                    {{1, 2, 6, 5}},
                                                    {{2, 3, 7, 6}},
                                                    {{3, 0, 4, 7}}}};

vtkNew<vtkActor> drawCube(double x, double y, double z, double h){
	std::array<std::array<double, 3>, 8> pts;
	pts[0] = {x,y,z};       //0,0,0
	pts[1] = {x+h, y, z };  //1,0,0
	pts[2] = {x+h, y+h, z}; //1,1,0
	pts[3] = {x, y+1, z};   //0,1,0
	pts[4] = {x, y, z+h};   //0,0,1
	pts[5] = {x+h,y,z+h};   //1,0,1
	pts[6] = {x+h,y+h,z+h}; //1,1,1
	pts[7] = {x,y+h,z+h};   //0,1,1

	vtkNew<vtkPolyData> cube;
	vtkNew<vtkPoints> points;
	vtkNew<vtkCellArray> polys;
	vtkNew<vtkFloatArray> scalars;

	for (auto i = 0ul; i < pts.size(); ++i) {
		points->InsertPoint(i, pts[i].data());
		scalars->InsertTuple1(i, i);
	}
	for (auto&& i : ordering) {
		polys->InsertNextCell(vtkIdType(i.size()), i.data());
	}

	// We now assign the pieces to the vtkPolyData.
	cube->SetPoints(points);
	cube->SetPolys(polys);
	cube->GetPointData()->SetScalars(scalars);

	// Now we'll look at it.
	vtkNew<vtkPolyDataMapper> cubeMapper;
	cubeMapper->SetInputData(cube);
	cubeMapper->SetScalarRange(cube->GetScalarRange());
	vtkNew<vtkActor> cubeActor;
	cubeActor->SetMapper(cubeMapper);
	cubeActor->GetProperty()->SetRepresentationToWireframe();

	return cubeActor;
}

vtkNew<vtkActor> drawPoint(double x, double y, double z) {
	vtkNew<vtkNamedColors> namedColors;

	vtkNew<vtkPoints> points;
	points->InsertNextPoint(x, y, z);

	vtkNew<vtkPolyData> pointsPolydata;
	pointsPolydata->SetPoints(points);


	vtkNew<vtkPolyData> polyData;

	polyData->SetPoints(points);

	vtkNew<vtkVertexGlyphFilter> glyphFilter;
	glyphFilter->SetInputData(polyData);
	glyphFilter->Update();
	

  
	// Visualization
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(glyphFilter->GetOutputPort());

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(10);
	actor->GetProperty()->SetColor(namedColors->GetColor3d("Olive").GetData());

	return actor;
}


int main(int argc, char **argv)
{
  vtkNew<vtkNamedColors> colors;

  vtkNew<vtkActor> cubeActor1 = drawCube(0, 0, 0, 1);
  vtkNew<vtkActor> cubeActor2 = drawCube(1, 1, 1, 1);

  // The usual rendering stuff.
  vtkNew<vtkCamera> camera;
  camera->SetPosition(1, 1, 1);
  camera->SetFocalPoint(0, 0, 0);

  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(renderer);
  renWin->SetWindowName("Cube");

  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);

  renderer->AddActor(cubeActor1);
  renderer->AddActor(cubeActor2);


  renderer->AddActor(drawPoint(0.5,0.5,0.5));
  renderer->SetActiveCamera(camera);
  renderer->ResetCamera();
  renderer->SetBackground(colors->GetColor3d("Cornsilk").GetData());

  renWin->SetSize(600, 600);

  // interact with data
  renWin->Render();
  iren->Start();

  return EXIT_SUCCESS;
}