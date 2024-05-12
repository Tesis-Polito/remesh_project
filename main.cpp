#include <igl/opengl/glfw/Viewer.h>
#include <igl/frame_field_deformer.h>
#include <igl/readOBJ.h>
#include <igl/avg_edge_length.h>
#include <iostream>
#include <fstream>
#include <igl/frame_to_cross_field.h>
#include <igl/find_cross_field_singularities.h>
#include <igl/per_face_normals.h>
#include <eigen3/Eigen/Core>
#include <igl/find_cross_field_singularities.h>



using namespace std;
using namespace Eigen;
using namespace igl;

void read_ff(string path, Eigen::MatrixXd& ff1, Eigen::MatrixXd& ff2);
void get_perp_vector(Eigen::MatrixXd& V1, Eigen::MatrixXd& N, Eigen::MatrixXd& P);


int main(int argc, char *argv[])
{
  //read the ground truth
  Eigen::MatrixXd V1;
  Eigen::MatrixXd bry;
  Eigen::MatrixXi F1;
  igl::readOBJ("/home/fabio/Documents/Automatic-field-remesher/objs/spheres_unclean.obj", V1, F1);

  
  // read presaved frame field
  Eigen::MatrixXd ff1;
  Eigen::MatrixXd ff2;
  read_ff("/home/fabio/Documents/testing_ff.ffpy",ff1,ff2);
  
  // frame field deformer (to get the singularities)
  Eigen::MatrixXd V1_d;
  Eigen::MatrixXd ff1_d;
  Eigen::MatrixXd ff2_d;
  Eigen::MatrixXd ff1_b;
  Eigen::MatrixXd ff2_b;

  Eigen::MatrixXd N;
  Eigen::VectorXi Singularities_id;
  Eigen::VectorXi isSingularity;




  frame_field_deformer(V1,F1,ff1, ff2, V1_d, ff1_d, ff2_d,100,0.5);

  per_face_normals(V1_d, F1, N);

  frame_to_cross_field(V1_d,F1,ff1_d,ff2_d, ff1_b); //ff1_b cross field vector, we need to rotate it 90`

  get_perp_vector(ff1_b, N,ff2_b);//we get the perpendicular vector to create the cross field (assuming normals with 1 norm)

  find_cross_field_singularities(V1_d, F1, ff1_b, ff2_b, isSingularity, Singularities_id, false);

 
  igl::barycenter(V1,F1,bry);
  // igl::barycenter(V1_d,F1,bry);

  double avg = avg_edge_length(V1,F1);
  // double avg = avg_edge_length(V1_d,F1);

  // // Plot the mesh
  const RowVector3d red(0.8,0.2,0.2),blue(0.2,0.2,0.8);
  igl::opengl::glfw::Viewer viewer;
  // viewer.data().set_mesh(V1_d, F1);
  viewer.data().set_mesh(V1, F1);

  viewer.data().set_face_based(true);
  // // //Cross field deform
  // viewer.data().add_edges(bry - ff1_b *avg /2, bry + ff1_b *avg/2, red);
  // viewer.data().add_edges(bry - ff2_b *avg /2, bry + ff2_b *avg/2, blue);

  // viewer.data().add_edges(bry - ff1_d *avg /2, bry + ff1_d *avg/2, red);
  // viewer.data().add_edges(bry - ff2_d *avg /2, bry + ff2_d *avg/2, blue);

  viewer.data().add_edges(bry - ff1 *avg /2, bry + ff1 *avg/2, red);
  viewer.data().add_edges(bry - ff2 *avg /2, bry + ff2 *avg/2, blue);

  for(unsigned i=0; i< V1_d.rows();i++){
    if (Singularities_id(i) < -0.001){
      // viewer.data().add_points(V1_d.row(i),RowVector3d(0,0,1));
      viewer.data().add_points(V1.row(i),RowVector3d(0,0,1));
    }

    else if (Singularities_id(i) > 0.001){
      // viewer.data().add_points(V1_d.row(i),RowVector3d(1,0,0));
      viewer.data().add_points(V1.row(i),RowVector3d(1,0,0));
    }
      

  }

  viewer.launch();
}

void read_ff(string path, Eigen::MatrixXd& ff1, Eigen::MatrixXd&ff2){
    std::ifstream fp (path);
    std::string line;
    double fx,fy,fz;
    double gx,gy,gz;
    Vector3d vec, vec1;
    int i=0;
    while (fp)
    {
      std::getline(fp, line);

      sscanf(line.c_str(), "[ %lf %lf %lf ] [ %lf %lf %lf ]", &fx,&fy,&fz, &gx, &gy, &gz);
      
      ff1.conservativeResize(ff1.rows()+1 , 3);
      ff1(i,0)=fx;
      ff1(i,1)=fy;
      ff1(i,2)=fz;
      ff2.conservativeResize(ff2.rows()+1 , 3);
      ff2(i)= gx;
      ff2(i,1)=gy;
      ff2(i,2)=gz;
      i++;
  
    }
    ff1.conservativeResize(ff1.rows()-1 , 3);
    ff2.conservativeResize(ff2.rows()-1 , 3);


    return;
}


void get_perp_vector(Eigen::MatrixXd& V1, Eigen::MatrixXd& N, Eigen::MatrixXd& P){

    Eigen::Vector3d vector1(0,0,0);
    Eigen::Vector3d normalv(0,0,0);

    P.conservativeResize(V1.rows(), 3);
    
    for (int i = 0; i < V1.rows(); i++)
    {
        vector1[0]= V1(i,0);
        vector1[1]=V1(i,1);
        vector1[2]=V1(i,2);

        normalv[0]=N(i,0);
        normalv[1]=N(i,1);
        normalv[2]=N(i,2);

        vector1= normalv.cross(vector1);

        P(i,0)=vector1[0];
        P(i,1)=vector1[1];
        P(i,2)=vector1[2];

    }

}