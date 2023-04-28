#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDataSetMapper.h>
#include <vtkImageData.h>
#include <vtkUnsignedCharArray.h>

#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

//Clase Punto
class Point {
public:
    double x, y, z;
    Point(double x_, double y_, double z_) {
        x = x_;
        y = y_;
        z = z_;
    }
    Point() {}

    Point operator-(const Point& otroPunto) {
        return Point(x - otroPunto.x, y - otroPunto.y, z - otroPunto.z);
    }
};

//Clase Distancia
class Distancia {
private:
    double dis;

public:
    Distancia(const Point& A, Point B) {
        dis += pow(A.x - B.x, 2.0);
        dis += pow(A.y - B.y, 2.0);
        dis += pow(A.z - B.z, 2.0);
        dis = sqrt(dis);
    }
    ~Distancia() {}

    double GetDis() {
        return dis;
    }
};

Point getMinPoint(Point*, int);
Point getMaxPoint(Point*, int);
double max(Point);

//Clase Octree
class Octree {
public:
    Octree* children[8] = { nullptr };
    Point* points;

    int K; //G granularidad
    int nPoints;

    Point leftBottom;
    double h;
    bool isLeaf = true; //si es hoja

    Octree() {}
    Octree(int K) {
        this->K = K;
        this->nPoints = 0;
        points = new Point[K];
    }
    ~Octree() {}

    // TDD: Test Driven Development
    // Ingresar el primer punto
    // Ingresar n puntos hasta antes que explote
    // Ingresar el punto k para que explote
    void insert(const Point& point) {
        if (isLeaf && nPoints == 0) {
            points[nPoints++] = point;
            h = 0;
            leftBottom = point;
        }
        else if (isLeaf && nPoints < K) {
            points[nPoints++] = point;
            Point minPoint = getMinPoint(points, nPoints);
            Point maxPoint = getMaxPoint(points, nPoints);
            h = max(minPoint - maxPoint);
            leftBottom = minPoint;
        }
        else {
            // buscar region que pertenece e insertarlo
            isLeaf = false;
            insertTypeOctree(point);
        }
    }

    void insertTypeOctree(const Point& point) {
        double midX = leftBottom.x + h / 2;
        double midY = leftBottom.y + h / 2;
        double midZ = leftBottom.z + h / 2;

        for (int i = 0; i < nPoints; i++) {
            if (points[i].x <= midX) {
                if (points[i].y <= midY) {
                    if (points[i].z <= midZ) {
                        if (!children[0])
                            children[0] = new Octree(K);

                        children[0]->insert(points[i]);
                    }
                    else {
                        if (!children[4])
                            children[4] = new Octree(K);

                        children[4]->insert(points[i]);
                    }
                }
                else {
                    if (points[i].z <= midZ) {
                        if (!children[2])
                            children[2] = new Octree(K);

                        children[2]->insert(points[i]);
                    }
                    else {
                        if (!children[6])
                            children[6] = new Octree(K);

                        children[6]->insert(points[i]);
                    }
                }
            }
            else {
                if (points[i].y <= midY) {
                    if (points[i].z <= midZ) {
                        if (!children[1])
                            children[1] = new Octree(K);

                        children[1]->insert(points[i]);
                    }
                    else {
                        if (!children[5])
                            children[5] = new Octree(K);

                        children[5]->insert(points[i]);
                    }
                }
                else {
                    if (points[i].z <= midZ) {
                        if (!children[3])
                            children[3] = new Octree(K);

                        children[3]->insert(points[i]);
                    }
                    else {
                        if (!children[7])
                            children[7] = new Octree(K);

                        children[7]->insert(points[i]);
                    }
                }
            }
        }

        points = nullptr;
        nPoints = 0;

        if (point.x <= midX) {
            if (point.y <= midY) {
                if (point.z <= midZ) {
                    if (!children[0])
                        children[0] = new Octree(K);

                    children[0]->insert(point);
                }
                else {
                    if (!children[4])
                        children[4] = new Octree(K);

                    children[4]->insert(point);
                }
            }
            else {
                if (point.z <= midZ) {
                    if (!children[2])
                        children[2] = new Octree(K);

                    children[2]->insert(point);
                }
                else {
                    if (!children[6])
                        children[6] = new Octree(K);

                    children[6]->insert(point);
                }
            }
        }
        else {
            if (point.y <= midY) {
                if (point.z <= midZ) {
                    if (!children[1])
                        children[1] = new Octree(K);

                    children[1]->insert(point);
                }
                else {
                    if (!children[5])
                        children[5] = new Octree(K);

                    children[5]->insert(point);
                }
            }
            else {
                if (point.z <= midZ) {
                    if (!children[3])
                        children[3] = new Octree(K);

                    children[3]->insert(point);
                }
                else {
                    if (!children[7])
                        children[7] = new Octree(K);

                    children[7]->insert(point);
                }
            }
        }
    }

    //caso basico: hoja
    //caso no es hoja
    bool exist(const Point& point) {
        if (isLeaf) {
            for (int i = 0; i < nPoints; i++) {
                if (points[i].x == point.x && points[i].y == point.y && points[i].z == point.z)
                    return true;
            }
            return false;
        }
        else {
            double midX = leftBottom.x + h / 2;
            double midY = leftBottom.y + h / 2;
            double midZ = leftBottom.z + h / 2;

            if (point.x <= midX) {
                if (point.y <= midY) {
                    if (point.z <= midZ)
                        children[0]->exist(point);
                    else
                        children[4]->exist(point);
                }
                else {
                    if (point.z <= midZ)
                        children[2]->exist(point);
                    else
                        children[6]->exist(point);
                }
            }
            else {
                if (point.y <= midY) {
                    if (point.z <= midZ)
                        children[1]->exist(point);
                    else
                        children[5]->exist(point);
                }
                else {
                    if (point.z <= midZ)
                        children[3]->exist(point);
                    else
                        children[7]->exist(point);
                }
            }
        }
    }

    Point closestPoint(const Point& point) {
        if (isLeaf) {
            double dist = (Distancia(point, points[0]).GetDis());
            Point closest = points[0];
            for (int i = 1; i < nPoints; i++) {
                Distancia(point, points[i]);
                double distTmp = (Distancia(point, points[0]).GetDis());
                if (distTmp < dist) {
                    dist = distTmp;
                    closest = points[i];
                }
            }
            return closest;
        }
        else {
            double midX = leftBottom.x + h / 2;
            double midY = leftBottom.y + h / 2;
            double midZ = leftBottom.z + h / 2;

            if (point.x <= midX) {
                if (point.y <= midY) {
                    if (point.z <= midZ)
                        children[0]->closestPoint(point);
                    else
                        children[4]->closestPoint(point);
                }
                else {
                    if (point.z <= midZ)
                        children[2]->closestPoint(point);
                    else
                        children[6]->closestPoint(point);
                }
            }
            else {
                if (point.y <= midY) {
                    if (point.z <= midZ)
                        children[1]->closestPoint(point);
                    else
                        children[5]->closestPoint(point);
                }
                else {
                    if (point.z <= midZ)
                        children[3]->closestPoint(point);
                    else
                        children[7]->closestPoint(point);
                }
            }
        }
    }
};

Point getMinPoint(Point* points, int nPoints) {
    Point minPoint = points[0];
    for (int i = 1; i < nPoints; i++) {
        if (points[i].x < minPoint.x)
            minPoint.x = points[i].x;
        if (points[i].y < minPoint.y)
            minPoint.y = points[i].y;
        if (points[i].z < minPoint.z)
            minPoint.z = points[i].z;
    }
    return minPoint;
}

Point getMaxPoint(Point* points, int nPoints) {
    Point minPoint = points[0];
    for (int i = 1; i < nPoints; i++) {
        if (points[i].x > minPoint.x)
            minPoint.x = points[i].x;
        if (points[i].y > minPoint.y)
            minPoint.y = points[i].y;
        if (points[i].z > minPoint.z)
            minPoint.z = points[i].z;
    }
    return minPoint;
}

double max(Point point) {
    if (abs(point.x) > abs(point.y) && abs(point.x) > abs(point.z))
        return abs(point.x);
    else if (abs(point.y) > abs(point.z))
        return abs(point.y);
    else return abs(point.z);
}


// The ordering of the corner points on each face.
std::array<std::array<vtkIdType, 4>, 6> ordering = {{{{0, 3, 2, 1}},
                                                    {{4, 5, 6, 7}},
                                                    {{0, 1, 5, 4}},
                                                    {{1, 2, 6, 5}},
                                                    {{2, 3, 7, 6}},
                                                    {{3, 0, 4, 7}}}};

vtkNew<vtkActor> drawCube(double x, double y, double z, double h) {
    std::array<std::array<double, 3>, 8> pts;
    pts[0] = {x,y,z};       //0,0,0
    pts[1] = {x+h, y, z };  //1,0,0
    pts[2] = {x+h, y+h, z}; //1,1,0
    pts[3] = {x,y+h,z};   //0,1,0
    pts[4] = {x,y,z+h};   //0,0,1
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

void cubosCoords(const Octree& octree, ofstream& archivo){
    if(octree.isLeaf){
         if(octree.h > -0.09 && octree.h < 0.09){
            archivo<<octree.leftBottom.x<<" "<<octree.leftBottom.y<<" "<< octree.leftBottom.z<<" "<<octree.h<<endl;
         }
    }
    else{
        for(int i = 0; i < 8; i++){
            if(octree.children[i] != nullptr){
                cubosCoords(*octree.children[i], archivo);
            }
        }
    }
}


int main() {
//Crear el Octree e insertar los puntos de dama_octal.txt
    Octree* root = new Octree(2);

    ifstream dama("dama_octal.txt");
    double x,y,z,h;
    while(dama>>x>>y>>z){
        root->insert(Point(x,y,z));
    }
    dama.close();

//Colocar x,y,z,h de las hojas en Octree.txt
    ofstream otCubos("Octree.txt");
    cubosCoords(*root, otCubos);
    otCubos.close();


//Leer cada cubo en Octree.txt y lo guarda en cubeActor
    vtkNew<vtkNamedColors> colors;
    vector<vtkNew<vtkActor>> cubeActor;

    ifstream cubos("Octree.txt");

    while(cubos>>x>>y>>z>>h){
        cubeActor.push_back(drawCube(x, y, z, h));
    }
    cubos.close();

    // The usual rendering stuff.
    vtkNew<vtkCamera> camera;
    camera->SetPosition(-5, -20, 15);
    camera->SetFocalPoint(0, 0, 0);

    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renWin;
    renWin->AddRenderer(renderer);
    renWin->SetWindowName("Octree");

    vtkNew<vtkRenderWindowInteractor> iren;
    iren->SetRenderWindow(renWin);

    for (int i = 0; i < cubeActor.size(); i++){
        renderer->AddActor(cubeActor[i]);
    }

    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();
    renderer->SetBackground(colors->GetColor3d("Cornsilk").GetData());

    renWin->SetSize(1800, 980);

    // interact with data
    renWin->Render();
    iren->Start();

    //###############################################
    cout << "Precione enter para continuar ... ";
    cin.get();

    return 0;
}