
#ifndef UNSTRUCTURED_MESH_H
#define UNSTRUCTURED_MESH_H
///////////////////////////////////////////////////////////////////////
//                  Unstructured_Mesh.h
//                  MCG4139 Winter 2021
///////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <Eigen/Core>
#include <vector>
#include <array>
#include <unordered_map>

///////////////////////////////////////////////////////////////////////
//                    Some useful types
///////////////////////////////////////////////////////////////////////

//structure to store centroid position
struct centroid_volume {
  double x;
  double y;
  double z;
};

//structure to store neighbour of each node
struct neighbour {
  int ne0=-1;
  int ne1=-1;
  int ne2=-1;
  int ne3=-1;
};

//structure to store properties of each node
struct volume {
  int node0;
  int node1;
  int node2;
  int node3;
  centroid_volume centroid;
  neighbour neighbours;
};

//structure to store properties of each triangle (face)
struct triangle {
  int node0;
  int node1;
  int node2;
    int v_up;
    int v_down;
    Eigen::Vector3d n_hat;
    double area;
    centroid_volume centroid;

};

using Node3D = Eigen::Vector3d;
using TriangleArray = std::vector<triangle>;


// Function to add neighbours to each node
void add_neighbours(volume &the_neighbour, volume &c, int neighbour_position, int current_position, int neighbour_neighbour_index, int current_neighbour_index) {
    if (current_neighbour_index == 1){
        c.neighbours.ne0 = neighbour_position;
    }else if (current_neighbour_index == 2){
        c.neighbours.ne1 = neighbour_position;
    }else if (current_neighbour_index == 3){
        c.neighbours.ne2 = neighbour_position;
    }else {
        c.neighbours.ne3 = neighbour_position;
    }
    if (neighbour_neighbour_index == 1){
        the_neighbour.neighbours.ne0 = current_position;
    }else if (neighbour_neighbour_index == 2){
        the_neighbour.neighbours.ne1 = current_position;
    }else if (neighbour_neighbour_index == 3){
        the_neighbour.neighbours.ne2 = current_position;
    }else {
        the_neighbour.neighbours.ne3 = current_position;
    }

}


////find centroid
auto compute_centroid(const Node3D& n0, const Node3D& n1, const Node3D& n2, const Node3D& n3) {
    centroid_volume centeroid;
    centeroid.x = (n0.x()+n1.x()+n2.x()+n3.x())/4;
    centeroid.y = (n0.y()+n1.y()+n2.y()+n3.y())/4;
    centeroid.z = (n0.z()+n1.z()+n2.z()+n3.z())/4;
  return centeroid;
}


// function to find the centroid of a traingle
auto compute_centroid_triangle(const Node3D& n0, const Node3D& n1, const Node3D& n2) {
    centroid_volume centeroid;
    centeroid.x = (n0.x()+n1.x()+n2.x())/3;
    centeroid.y = (n0.y()+n1.y()+n2.y())/3;
    centeroid.z = (n0.z()+n1.z()+n2.z())/3;
  return centeroid;
}

// calculating the n_hat and area for each triangle
auto n_hat_and_area(const Node3D& n0, const Node3D& n1, const Node3D& n2) {
    const Node3D v1 = n1-n0;
    const Node3D v2 = n2-n0;
    auto area = 0.5*sqrt((v1.y()*v2.z()-v1.z()*v2.y())*(v1.y()*v2.z()-v1.z()*v2.y())+(v1.z()*v2.x()-v1.x()*v2.z())*(v1.z()*v2.x()-v1.x()*v2.z())+(v1.x()*v2.y()-v1.y()*v2.x())*(v1.x()*v2.y()-v1.y()*v2.x()));
  Node3D n_hat;

    n_hat.x() = 0.5*(v1.y()*v2.z() - v1.z()*v2.y())/area;
    n_hat.y() = 0.5*(v1.z()*v2.x() - v1.x()*v2.z())/area;
    n_hat.z() = 0.5*(v1.x()*v2.y() - v1.y()*v2.x())/area;

    return std::make_pair(n_hat, area);
}

///////////////////////////////////////////////////////////////////////
//                    Read gmsh file
///////////////////////////////////////////////////////////////////////
auto read_gmsh_file(const std::string& filename) {

    

  std::ifstream fin(filename);
  if(!fin) {
    throw std::runtime_error("Cannot open file: " + filename + ".");
  }

  ///////////////////////////////////////////////////////
  // Lambda function to consume lines that are expected
  // but should be ignored
  auto expect_line = [filename] (std::ifstream& fin, const std::string& expected) {
    std::string s;
    do {
      std::getline(fin,s);
    } while (s.empty());

    if(s != expected) {
      throw std::runtime_error("Error reading file: " + filename + ".\n" +
                               "Expected \"" + expected + "\", but got \"" + s +
                               "\".");
    }
  };


  std::cout << "Reading gmsh file: " << filename << '\n';

  ///////////////////////////////////////////////////////
  // Read file
  expect_line(fin, "$MeshFormat");
  expect_line(fin, "2.2 0 8");
  expect_line(fin, "$EndMeshFormat");
  expect_line(fin, "$Nodes");

  int number_of_nodes;
  fin >> number_of_nodes;

  // reading nodes indexes
  std::vector<Node3D> nodes(number_of_nodes);
  for(int i = 0; i < number_of_nodes; ++i) {
    int    dummy_index;
    fin >> dummy_index
        >> nodes[i].x()
        >> nodes[i].y()
        >> nodes[i].z();
    if(dummy_index != i+1) {
      throw std::runtime_error("Error with node index.");
    }
  }

  expect_line(fin, "$EndNodes");
  expect_line(fin, "$Elements");

  int number_of_elements; //not all will be cells
  fin >> number_of_elements;

  std::vector<volume> volumes;
  std::vector<triangle> triangles;
  std::unordered_map< std::string ,std::array<int,3>> mymap;
    
  int j =-1;
    
  for(int i = 0; i < number_of_elements; ++i) {
    std::string s;
    do {
      std::getline(fin,s);
    } while (s.empty());

    std::istringstream ss(s);

    int element_num;
    ss >> element_num;
    if(element_num - 1 != i) {
      throw std::runtime_error("Error reading element number.");
    }

    int element_type;
    ss >> element_type;
     
    if(element_type == 4) {
      //tetra cell
      j = j+1;
      int dummy;
      int n[4];
      volume c;

      ss >> dummy >> dummy >> dummy >> n[0] >> n[1] >> n[2] >> n[3];
      std::sort(n, n+4);
        
      c.node0 = n[0]-1;
      c.node1 = n[1]-1;
      c.node2 = n[2]-1;
      c.node3 = n[3]-1;


      centroid_volume c2;
      c2 = compute_centroid(nodes[c.node0], nodes[c.node1], nodes[c.node2], nodes[c.node3]);
      c.centroid = c2;
        int neighbour_position;
        int neighbour_neighbour_index;
   
        // storing all the nodes in an unordered map
        // for each triangle, find its neighbours
        std::string key = std::to_string(c.node0) +","+std::to_string(c.node1)+","+std::to_string(c.node2);
        if (mymap.find(key) == mymap.end()){
            
            triangle c3;
            c3.node0=c.node0;
            c3.node1=c.node1;
            c3.node2=c.node2;
            triangles.push_back(c3);
            
            mymap[key] = {j,1,-1};
        }else {
            mymap[key][2] = j;
            neighbour_position = mymap[key][0];
            neighbour_neighbour_index = mymap[key][1];
            add_neighbours(volumes[neighbour_position], c, neighbour_position, j,neighbour_neighbour_index, 1);
        }
        key = std::to_string(c.node1) +","+std::to_string(c.node2)+","+std::to_string(c.node3);
        if (mymap.find(key) == mymap.end()){
            
            triangle c3;
            c3.node0=c.node1;
            c3.node1=c.node2;
            c3.node2=c.node3;
            triangles.push_back(c3);
            
            mymap[key] = {j,2,-1};
        }else {
            mymap[key][2] = j;
            neighbour_position = mymap[key][0];
            neighbour_neighbour_index = mymap[key][1];
            add_neighbours(volumes[neighbour_position], c, neighbour_position, j,neighbour_neighbour_index, 2);
        }
        key = std::to_string(c.node0) +","+std::to_string(c.node2)+","+std::to_string(c.node3);
        if (mymap.find(key) == mymap.end()){
            
            triangle c3;
            c3.node0=c.node0;
            c3.node1=c.node2;
            c3.node2=c.node3;
            triangles.push_back(c3);
            
            mymap[key] = {j,3,-1};
        }else {
            mymap[key][2] = j;
            neighbour_position = mymap[key][0];
            neighbour_neighbour_index = mymap[key][1];
            add_neighbours(volumes[neighbour_position], c, neighbour_position, j,neighbour_neighbour_index, 3);
        }
        key = std::to_string(c.node0) +","+std::to_string(c.node1)+","+std::to_string(c.node3);
        if (mymap.find(key) == mymap.end()){
            
            triangle c3;
            c3.node0=c.node0;
            c3.node1=c.node1;
            c3.node2=c.node3;
            triangles.push_back(c3);
            
            mymap[key] = {j,4,-1};
        }else {
            mymap[key][2] = j;
            neighbour_position = mymap[key][0];
            neighbour_neighbour_index = mymap[key][1];
            add_neighbours(volumes[neighbour_position], c, neighbour_position, j,neighbour_neighbour_index, 4);
        }
      volumes.push_back(c);
    }

  }


  expect_line(fin, "$EndElements");

  std::cout << "done.\n";
    
    // for each triangle calculate its properites, including n_hat, area, and centroid
    for(auto& triangle1 : triangles) {
        double area;
        Node3D n_hat;
        std::tie(n_hat,area) = n_hat_and_area(nodes[triangle1.node0], nodes[triangle1.node1], nodes[triangle1.node2]);
        triangle1.n_hat = n_hat;
        triangle1.area = area;
        std::string key = std::to_string(triangle1.node0) +","+std::to_string(triangle1.node1)+","+std::to_string(triangle1.node2);
        int v_up_index = mymap[key][0];
        int v_down_index = mymap[key][2];
        int fourth_node = mymap[key][1];
        volume v_up = volumes[v_up_index];
        if (fourth_node == 1){
            fourth_node = v_up.node3;
        } else if (fourth_node == 2){
            fourth_node = v_up.node0;
        } else if (fourth_node == 3){
            fourth_node = v_up.node1;
        } else if (fourth_node == 4){
            fourth_node = v_up.node2;
        }
        Node3D check_if_same_direction = nodes[fourth_node] - nodes[triangle1.node0];
        
        if (check_if_same_direction.x()*n_hat.x() + check_if_same_direction.y()*n_hat.y()+ check_if_same_direction.z()*n_hat.z() >=0 ){

            triangle1.v_up = v_up_index;
            triangle1.v_down = v_down_index;
        }else {

            triangle1.v_down = v_up_index;
            triangle1.v_up = v_down_index;
        }
        
        centroid_volume c2;
        c2 = compute_centroid_triangle(nodes[triangle1.node0], nodes[triangle1.node1], nodes[triangle1.node2]);
        triangle1.centroid = c2;
    }
   
  return std::make_tuple(nodes, volumes, triangles);
}


#endif //#ifndef UNSTRUCTURED_MESH_H


