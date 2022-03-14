#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <iostream>
#include <filesystem>
#include "scene.h"
#include "line_cylinders.h"


Eigen::MatrixXd V;
Eigen::MatrixXi F;
igl::opengl::glfw::Viewer mgpViewer;

float currTime = 0;
bool animationHack;  //fixing the weird camera bug in libigl
//initial values
float timeStep = 0.02;
float CRCoeff= 1.0;


double tolerance = 10e-3;
int maxIterations=10000;

Scene scene;

Eigen::MatrixXd platV;
Eigen::MatrixXi platF;
Eigen::MatrixXi platT;
Eigen::RowVector3d platCOM;
Eigen::RowVector4d platOrientation;

void createPlatform()
{
  double platWidth=100.0;
  platCOM<<0.0,-5.0,-0.0;
  platV.resize(9,3);
  platF.resize(12,3);
  platT.resize(12,4);
  platV<<-platWidth,0.0,-platWidth,
  -platWidth,0.0,platWidth,
  platWidth,0.0,platWidth,
  platWidth,0.0, -platWidth,
  -platWidth,-platWidth/10.0,-platWidth,
  -platWidth,-platWidth/10.0,platWidth,
  platWidth,-platWidth/10.0,platWidth,
  platWidth,-platWidth/10.0, -platWidth,
  0.0,-platWidth/20.0, 0.0;
  platF<<0,1,2,
  2,3,0,
  6,5,4,
  4,7,6,
  1,0,5,
  0,4,5,
  2,1,6,
  1,5,6,
  3,2,7,
  2,6,7,
  0,3,4,
  3,7,4;
  
  platOrientation<<1.0,0.0,0.0,0.0;
  
  platT<<platF, VectorXi::Constant(12,8);
  
  
}

void updateMeshes(igl::opengl::glfw::Viewer &viewer)
{
  RowVector3d platColor; platColor<<0.8,0.8,0.8;
  RowVector3d meshColor; meshColor<<0.8,0.2,0.2;
  viewer.core().align_camera_center(scene.meshes[0].currV);
  for (int i=0;i<scene.meshes.size();i++){
    viewer.data_list[i].clear();
    viewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
    viewer.data_list[i].set_face_based(true);
    viewer.data_list[i].set_colors(meshColor);
    viewer.data_list[i].show_lines=false;
  }
  viewer.data_list[0].show_lines=false;
  viewer.data_list[0].set_colors(platColor.replicate(scene.meshes[0].F.rows(),1));
  viewer.data_list[0].set_face_based(true);
  //viewer.core.align_camera_center(scene.meshes[0].currV);
  
  //updating constraint viewing
  MatrixXi constF;
  MatrixXd constV, constC;
  
  MatrixXd P1(scene.constraints.size(),3);
  MatrixXd P2(scene.constraints.size(),3);
  for (int i=0;i<scene.constraints.size();i++){
    P1.row(i)=scene.meshes[scene.constraints[i].m1].currV.row(scene.constraints[i].v1);
    P2.row(i)=scene.meshes[scene.constraints[i].m2].currV.row(scene.constraints[i].v2);
  }
  
  MatrixXd cyndColors=RowVector3d(1.0,1.0,0.0).replicate(P1.size(),1);
  
  double radius = 0.5;
  hedra::line_cylinders(P1,P2,radius,cyndColors,8,constV,constF,constC);
  viewer.data_list[scene.meshes.size()].set_mesh(constV, constF);
  viewer.data_list[scene.meshes.size()].set_face_based(true);
  viewer.data_list[scene.meshes.size()].set_colors(constC);
  viewer.data_list[scene.meshes.size()].show_lines=false;
  
  
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  if (key == ' ')
  {
    viewer.core().is_animating = !viewer.core().is_animating;
    if (viewer.core().is_animating)
      cout<<"Simulation running"<<endl;
    else
      cout<<"Simulation paused"<<endl;
    return true;
  }
  
  if (key == 'S')
  {
    if (!viewer.core().is_animating){
      scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations);
      currTime+=timeStep;
      updateMeshes(viewer);
      std::cout <<"currTime: "<<currTime<<std::endl;
      return true;
    }
  }
  return false;
}




bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
  using namespace Eigen;
  using namespace std;
  
  if (viewer.core().is_animating){
    if (!animationHack)
      scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations);
    else
      viewer.core().is_animating=false;
    animationHack=false;
    currTime+=timeStep;
    //cout <<"currTime: "<<currTime<<endl;
    updateMeshes(viewer);
  }
 
  
  return false;
}

class CustomMenu : public igl::opengl::glfw::imgui::ImGuiMenu
{
  
  virtual void draw_viewer_menu() override
  {
    // Draw parent menu
    //ImGuiMenu::draw_viewer_menu();
    
    // Add new group
    if (ImGui::CollapsingHeader("Algorithm Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
      ImGui::InputFloat("CR Coeff",&CRCoeff,0,0,"%.2f");
      
      
      if (ImGui::InputFloat("Time Step", &timeStep)) {
        mgpViewer.core().animation_max_fps = (((int)1.0/timeStep));
      }
    }
  }
};

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  const std::string DATA_DIR_STR = TOSTRING(DATA_DIR_M);
  std::vector<std::string> scene_files;

  for (const auto& file : std::filesystem::directory_iterator(DATA_DIR_STR))
  {
      if (!file.is_regular_file()) continue;
      if (!(file.path().extension().string() == ".txt")) continue;
      scene_files.push_back(file.path().filename().string());
  }

  cout << "\nScene files found in data folder:" << endl;
  for (int i = 0; i < scene_files.size(); i++)
  {
      cout << i << ") " << scene_files[i] << endl;
  }
  cout << endl;

  int sceneIndex = -1;
  while (sceneIndex < 0 || sceneIndex >= scene_files.size())
  {
      cout << "Please choose a scene index(0 - " << (scene_files.size() - 1) << ")" << endl;
      cin >> sceneIndex;
  }

  int constraintIndex = -1;
  while (constraintIndex < 0 || constraintIndex >= scene_files.size())
  {
      cout << "Please choose a constraint index(0 - " << (scene_files.size() - 1) << ")" << endl;
      cin >> constraintIndex;
  }

  //create platform
  createPlatform();
  scene.addMesh(platV, platF, platT, 10000.0, true, platCOM, platOrientation);
  
  //load scene from file
  scene.loadScene(DATA_DIR_STR, scene_files[sceneIndex], scene_files[constraintIndex]);

  scene.updateScene(0.0, CRCoeff, tolerance, maxIterations);
  
  // Viewer Settings
  for (int i=0;i<scene.meshes.size();i++){
    if (i!=0)
      mgpViewer.append_mesh();
    //mgpViewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
  }
  //mgpViewer.core.align_camera_center(scene.meshes[0].currV);
  
  //constraints mesh (for lines)
  mgpViewer.append_mesh();
  mgpViewer.callback_pre_draw = &pre_draw;
  mgpViewer.callback_key_down = &key_down;
  mgpViewer.core().is_animating = true;
  animationHack = true;
  mgpViewer.core().animation_max_fps = 50.;
  
  CustomMenu menu;
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  mgpViewer.plugins.push_back(&plugin);
  plugin.widgets.push_back(&menu);
  
  cout<<"Press [space] to toggle continuous simulation" << endl;
  cout<<"Press 'S' to advance time step-by-step"<<endl;
  
  updateMeshes(mgpViewer);
  mgpViewer.launch();
 
}
