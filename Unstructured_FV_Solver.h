#ifndef UNSTRUCTURED_FV_SOLVER_H
#define UNSTRUCTURED_FV_SOLVER_H
///////////////////////////////////////////////////////////////////////
//                  Unstructured_FV_Solver.h
///////////////////////////////////////////////////////////////////////
#include <cmath>
#include <string>
#include <functional>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "Unstructured_Mesh.h"
#include <array>
#include <vector>
#include <math.h>


//structure to store temperature slopes
struct ABC {
  double A;
  double B;
  double C;
};


//invert a matrix
auto invert_matrix (Eigen::Matrix<double,3,3> &m)
{
    double det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
                 m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
                 m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    //det = 0, then no inverse
    if (fabs(det) <= 1e-8){
        m.fill(0.0);
        return m;
    }
    double invdet = 1 / det;

    Eigen::Matrix<double,3,3> minv; // inverse of matrix m
    minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
    minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
    minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
    minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
    minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
    minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
    minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
    minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
    minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;
    return minv;
}

//adding information to the right hand side of the equation to solve for slopes
void add_to_RHS (volume &current_volume, double &current_T, volume &n_volume, double &n_T, Eigen::VectorXd &RHS){
    RHS(0) += 2.0*(n_T-current_T)*(n_volume.centroid.x-current_volume.centroid.x);
    RHS(1) += 2.0*(n_T-current_T)*(n_volume.centroid.y-current_volume.centroid.y);
    RHS(2) += 2.0*(n_T-current_T)*(n_volume.centroid.z-current_volume.centroid.z);
}
//adding information to the left hand side of the equation to solve for slopes
void add_to_matrix (volume &current_volume, volume &n_volume, Eigen::Matrix<double,3,3> &eq)
{
    eq(0,0) += 2.0*(n_volume.centroid.x-current_volume.centroid.x)*(n_volume.centroid.x-current_volume.centroid.x);
    eq(0,1) += 2.0*(n_volume.centroid.y-current_volume.centroid.y)*(n_volume.centroid.x-current_volume.centroid.x);
    eq(0,2) += 2.0*(n_volume.centroid.z-current_volume.centroid.z)*(n_volume.centroid.x-current_volume.centroid.x);
    eq(1,0) += 2.0*(n_volume.centroid.x-current_volume.centroid.x)*(n_volume.centroid.y-current_volume.centroid.y);
    eq(1,1) += 2.0*(n_volume.centroid.y-current_volume.centroid.y)*(n_volume.centroid.y-current_volume.centroid.y);
    eq(1,2) += 2.0*(n_volume.centroid.z-current_volume.centroid.z)*(n_volume.centroid.y-current_volume.centroid.y);
    eq(2,0) += 2.0*(n_volume.centroid.x-current_volume.centroid.x)*(n_volume.centroid.z-current_volume.centroid.z);
    eq(2,1) += 2.0*(n_volume.centroid.y-current_volume.centroid.y)*(n_volume.centroid.z-current_volume.centroid.z);
    eq(2,2) += 2.0*(n_volume.centroid.z-current_volume.centroid.z)*(n_volume.centroid.z-current_volume.centroid.z);
}

// find the inverse of the left hand side of the equation for each volume
void calc_eq (std::vector<volume> &volumes, std::vector<Eigen::Matrix<double,3,3>> &eq){
    Eigen::Matrix<double,3,3> inv_eq;
    for (int i = 0; i < volumes.size() ; ++i){
        inv_eq.fill(0.0);
        if (volumes[i].neighbours.ne0 != -1){
            add_to_matrix(volumes[i], volumes[volumes[i].neighbours.ne0], inv_eq );
        }
        if (volumes[i].neighbours.ne1 != -1){
            add_to_matrix(volumes[i],volumes[volumes[i].neighbours.ne1], inv_eq );
        }
        if (volumes[i].neighbours.ne2 != -1){
            add_to_matrix(volumes[i],volumes[volumes[i].neighbours.ne2], inv_eq );
        }
        if (volumes[i].neighbours.ne3 != -1){
            add_to_matrix(volumes[i],volumes[volumes[i].neighbours.ne3], inv_eq );
        }
        inv_eq = invert_matrix(inv_eq);
        eq.push_back(inv_eq);
    }
}

//compute the slopes for each volume by nultiplying the inverse of the left hand side matrix by the right hand side vector
auto comupte_A_B_C(std::vector<volume> &volumes, std::vector<double> &U,std::vector<Eigen::Matrix<double,3,3>> &eq){
    std::vector<ABC> full_solution;
    auto RHS = Eigen::VectorXd(3);

    for (int i = 0; i < volumes.size() ; ++i){
        RHS.fill(0.0);
        if (volumes[i].neighbours.ne0 != -1){
            add_to_RHS(volumes[i], U[i],volumes[volumes[i].neighbours.ne0], U[volumes[i].neighbours.ne0], RHS );
        }
        if (volumes[i].neighbours.ne1 != -1){
            add_to_RHS(volumes[i], U[i],volumes[volumes[i].neighbours.ne1], U[volumes[i].neighbours.ne1], RHS );
        }
        if (volumes[i].neighbours.ne2 != -1){
            add_to_RHS(volumes[i], U[i],volumes[volumes[i].neighbours.ne2], U[volumes[i].neighbours.ne2], RHS );
        }
        if (volumes[i].neighbours.ne3 != -1){
            add_to_RHS(volumes[i], U[i],volumes[volumes[i].neighbours.ne3], U[volumes[i].neighbours.ne3], RHS );
        }

        auto solution = Eigen::VectorXd(3);
        solution = eq[i]*RHS;
        ABC full;
        
        full.A = solution(0);
        full.B = solution(1);
        full.C = solution(2);
        
        full_solution.push_back(full);
    }
    return full_solution;
}

///////////////////////////////////////////////////////////////////////
//             Unstructured Finite-Volume Scheme
///////////////////////////////////////////////////////////////////////
class Unstructured_FV_Solver {
public:
  ///////////////////////////////////////////////////////////////////////
  //  Default Constructors, etc.
  Unstructured_FV_Solver() = default;
  Unstructured_FV_Solver(const Unstructured_FV_Solver&) = default;
  Unstructured_FV_Solver(Unstructured_FV_Solver&&) = default;
  Unstructured_FV_Solver& operator=(const Unstructured_FV_Solver&) = default;
  Unstructured_FV_Solver& operator=(Unstructured_FV_Solver&&) = default;

  ///////////////////////////////////////////////////////////////////////
  //  Constructor taking a filename and initial condition
  Unstructured_FV_Solver(const std::string& mesh_filename);

  ///////////////////////////////////////////////////////////////////////
  //  Solution vector in cell "i"
  auto& U(int i) {
    return Global_U[i];
  }

  ///////////////////////////////////////////////////////////////////////
  //  dUdt in cell "i"
  auto& dUdt(int i) {
    return Global_dUdt[i];
  }

  ///////////////////////////////////////////////////////////////////////
  //  number of cells, nodes, and edges.
  auto number_of_cells() {return volumes.size();}
  auto number_of_nodes() {return nodes.size();}

  ///////////////////////////////////////////////////////////////////////
  //  time march to time
  void time_march_to_time(double final_time, double CFL, std::vector<Eigen::Matrix<double,3,3>> &eq);

  ///////////////////////////////////////////////////////////////////////
  //  Make movie
  void make_movie(double final_time, double CFL, int num_frames, const std::string& filename_base);

  ///////////////////////////////////////////////////////////////////////
  //  write to VTK
  void write_to_vtk(const std::string& filename);

private:

  ///////////////////////////////////////////////////////////////////////
  //  Member variables
  std::vector<Node3D>                                  nodes;
  std::vector<volume>                                  volumes;
  std::vector<triangle>                                triangles;
  std::vector<double>                                  volume_of_volumes;
  std::vector<double>                                  Global_U;
  std::vector<double>                                  Global_dUdt;
  double                                               time;

    
   
  ///////////////////////////////////////////////////////////////////////
    //dot product of 2 vectors
    auto dot_product(const Node3D &v1, const Node3D &v2){
        return (v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z());
    }

    //cross product of 2 vectors
    auto cross_product(const Node3D &v1, const Node3D &v2){
        Node3D cross;
        cross.x() = v1.y()*v2.z() - v1.z()*v2.y();
        cross.y() = v1.z()*v2.x() - v1.x()*v2.z();
        cross.z() = v1.x()*v2.y() - v1.y()*v2.x();
        return cross;
    }
    //calculating the volume
    auto calc_volume(){
        volume_of_volumes.resize(number_of_cells());
        
        for(int i = 0; i < number_of_cells(); ++i) {
          const volume& cell = volumes[i];
            const Node3D v1 = nodes[cell.node1]-nodes[cell.node0];
            const Node3D v2 = nodes[cell.node2]-nodes[cell.node0];
            const Node3D v3 = nodes[cell.node3]-nodes[cell.node0];
            const Node3D cross = cross_product(v1, v2);
            double v = dot_product(cross, v3);
            volume_of_volumes[i] = fabs((1.0/6.0)*v);
        }
    }

};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//  Constructor

Unstructured_FV_Solver::Unstructured_FV_Solver(const std::string& mesh_filename) {

  std::tie(nodes,volumes, triangles) = read_gmsh_file(mesh_filename);

  Global_U.resize(number_of_cells());
  Global_dUdt.resize(number_of_cells());
    
    //initial condition
    for(int i = 0; i < number_of_cells(); ++i) {
        U(i) =20;
    }
    
    calc_volume();
  time = 0.0;
}




///////////////////////////////////////////////////////////////////////
//  time march to time

void Unstructured_FV_Solver::time_march_to_time(double final_time, double CFL, std::vector<Eigen::Matrix<double,3,3>> &eq) {

  constexpr double tolerance = 1.0e-12;
    double k =0.9; // diffusivity
    double Ta = 20; //ambiant temperature
    double h = 0.03; //convection coefficient
    double source = 100; //source in watts

  while(time < final_time - tolerance) {

      double dt = 2.5e-6;
      std::vector<ABC> full_solution = comupte_A_B_C(volumes, Global_U, eq);
    
      for(auto& entry : Global_dUdt) {
          entry = 0.0;
      }
      for(int i = 0; i < number_of_cells(); ++i) {
        dt = std::min(dt, fabs(0.5*CFL*pow(volume_of_volumes[i], 1.0/3.0))/k);
      }
      double dU_down = 0.0;
      double dU_up = 0.0;
      double U_down = 0.0;
      double U_up = 0.0;
      for(const auto& triangle1 : triangles) {
          
          
          double F;
          if(triangle1.v_down >= 0 && triangle1.v_up >= 0) {
            dU_down = full_solution[triangle1.v_down].A*triangle1.n_hat.x()+full_solution[triangle1.v_down].B*triangle1.n_hat.y()+full_solution[triangle1.v_down].C*triangle1.n_hat.z();
              
            dU_up = full_solution[triangle1.v_up].A*triangle1.n_hat.x()+full_solution[triangle1.v_up].B*triangle1.n_hat.y()+full_solution[triangle1.v_up].C*triangle1.n_hat.z();

              
              F = -0.5*k*dt*(dU_down+dU_up)-sqrt((k*dt)/M_PI)*(U(triangle1.v_up)-U(triangle1.v_down));
             
                  dUdt(triangle1.v_down) += -F*triangle1.area/volume_of_volumes[triangle1.v_down];
            
                  dUdt(triangle1.v_up) += F*triangle1.area/volume_of_volumes[triangle1.v_up];
              
              
          }else if(triangle1.v_down == -1) {
              
              U_up = U(triangle1.v_up) + full_solution[triangle1.v_up].A*(triangle1.centroid.x-volumes[triangle1.v_up].centroid.x)+full_solution[triangle1.v_up].B*(triangle1.centroid.y-volumes[triangle1.v_up].centroid.y)+full_solution[triangle1.v_up].C*(triangle1.centroid.z-volumes[triangle1.v_up].centroid.z);
              
              if(nodes[triangle1.node0].z() == 0 && nodes[triangle1.node1].z() == 0 && nodes[triangle1.node2].z() == 0){
                  F = source*dt;
              }else {
                  F = -h*(U_up-Ta)*dt;
              }
              dUdt(triangle1.v_up) += F*triangle1.area/volume_of_volumes[triangle1.v_up];
          }else if(triangle1.v_up == -1) {
              
              U_down = U(triangle1.v_down) + full_solution[triangle1.v_down].A*(triangle1.centroid.x-volumes[triangle1.v_down].centroid.x)+full_solution[triangle1.v_down].B*(triangle1.centroid.y-volumes[triangle1.v_down].centroid.y)+full_solution[triangle1.v_down].C*(triangle1.centroid.z-volumes[triangle1.v_down].centroid.z);
              if(nodes[triangle1.node0].z() == 0 && nodes[triangle1.node1].z() == 0 && nodes[triangle1.node2].z() == 0){
                  F = source*dt;
              }else {
                  F = -h*(U_down-Ta)*dt;
              }
              dUdt(triangle1.v_down) += F*triangle1.area/volume_of_volumes[triangle1.v_down];
            
            }
 
      }
      for(int i = 0; i < number_of_cells(); ++i) {
          
          U(i) += dUdt(i);
    
      }
      

    time += dt;

    std::cout << "Time = " << time
               << "    dt = " << dt << '\n';
  }

}


///////////////////////////////////////////////////////////////////////
//  Make movie
void Unstructured_FV_Solver::make_movie(double final_time,
                                                  double CFL,
                                                  int num_frames,
                                                  const std::string& filename_base) {

  const auto dt = final_time/static_cast<double>(num_frames-1);

  auto get_name = [&filename_base] (int i) {
    std::stringstream ss;
    ss << filename_base << std::setfill('0') << std::setw(10) << i << ".vtk";
    return ss.str();
  };
  std::vector<Eigen::Matrix<double,3,3>> eq;
  calc_eq(volumes, eq);
  write_to_vtk(get_name(0));

  for(int i = 1; i < num_frames; ++i) {
    std::cout << "Time Marching to time =" << dt*static_cast<double>(i) << '\n';
    time_march_to_time(dt*static_cast<double>(i), CFL, eq);
    write_to_vtk(get_name(i));
  }

}



///////////////////////////////////////////////////////////////////////
//  write to VTK

void Unstructured_FV_Solver::write_to_vtk(const std::string& filename) {

  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("error opening vtk file.");
  }

  fout << "# vtk DataFile Version 2.0\n"
       << "Unstructured Solver\n"
       << "ASCII\n"
       << "DATASET UNSTRUCTURED_GRID\n"
       << "POINTS " << number_of_nodes() << " double\n";

  for(const auto& node : nodes) {
    fout << node.x() << " " << node.y() << " " << node.z() << '\n';
  }

  fout << "\nCELLS " << number_of_cells() << " " << 5*number_of_cells() << '\n';
  for(const auto& cell : volumes) {
    fout << "4 " << cell.node0 << " " << cell.node1 << " " << cell.node2 << " " <<cell.node3 << '\n';
  }
////change number  "4"
  fout << "\nCELL_TYPES " << number_of_cells() << '\n';
  for(int i = 0; i < number_of_cells(); ++i) {
    fout << "10\n";
  }

  fout << "\nCELL_DATA " << number_of_cells()  << '\n'
       << "SCALARS variable1 double 1\n"
       << "LOOKUP_TABLE default\n";
  for(int i = 0; i < number_of_cells(); ++i) {
    fout << U(i) << '\n';
  }

}

#endif //#ifndef UNSTRUCTURED_FV_SOLVER_H
