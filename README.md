
# Unstructured Finite Volume Solver

This project implements a simple unstructured Finite Volume solver designed to simulate fluid dynamics or other physical processes on unstructured meshes. The solver reads a mesh file, runs a simulation, and generates an output movie to visualize the results.

# Project Overview

The solver operates on unstructured mesh grids, using a finite volume method to solve for physical quantities over the mesh. It is implemented in C++ and consists of the following key components:

1. Unstructured_Mesh: Manages the mesh data and operations related to the mesh structure.
   
2. Unstructured_FV_Solver: Implements the finite volume solver, handling the simulation logic and time-stepping.
   
3. Main Program (main.cpp): The entry point that creates the solver and starts the simulation.

# Classes and Their Properties

1. Unstructured_Mesh
The Unstructured_Mesh class is responsible for reading, storing, and managing the mesh data. It handles the mesh's topology and structure.

Key Properties:

	Nodes: Stores the list of nodes (vertices) in the mesh.
	Cells: Represents the cells (elements) of the mesh, which could be triangles, quadrilaterals, or other polygons depending on the mesh type.
	Faces: Stores information about the faces of the cells, which are important for defining boundary conditions.

Key Functions:

	read_gmsh_file: Reads the mesh data from a file.
	n_hat_and_area: calculate the n_hat and area for each cell.
	compute_centroid: calculate the centroid for each cell.
	add_neighbours: Add the neighboring cells of a given cell to the cell structure.
	This class ensures that the mesh data is properly structured and accessible for the solver to work on.

2. Unstructured_FV_Solver
The Unstructured_FV_Solver class is where the core functionality of the finite volume solver resides. It is responsible for setting up and running the simulation, updating the solution at each time step, and generating the output movie.

Key Properties:

	Mesh: An instance of the Unstructured_Mesh class. This provides the solver with the mesh data on which it operates.
	Time Step (dt): The time step used for numerical integration in the simulation.
	Simulation Parameters: Various parameters such as the total simulation time, number of frames to generate, and output file names are defined here.

Key Functions:

	time_march_to_time: Performs the time marching to the specified final_time, using the Courant-Friedrichs-Lewy (CFL) condition.
	make_movie: This function controls the creation of the output movie. It takes parameters such as the time step, frame interval, and the name of the output movie file. It runs the simulation, updates the solution, and generates images for each frame, which are then compiled into a movie.
	write_to_vtk: Writes the solution data to a VTK file for visualization in tools like ParaView.
	dot_product: Computes the dot product of two 3D vectors.
	cross_product: Computes the cross product of two 3D vectors.
	calc_volume(): Computes the volume for each mesh cell using the determinant of a 3x3 matrix formed by the vectors from the nodes of the volume.

3. Main Program (main.cpp)
The main program serves as the entry point for the project. It initializes the solver with a mesh file, calls the solver to run the simulation, and generates an animation.

Key Steps:

	Create Solver Instance: An instance of the Unstructured_FV_Solver class is created, with the mesh file (nofins.msh) passed to the constructor.
	Call make_movie: The make_movie function is called to run the simulation and generate the movie. It takes parameters like time step, frame interval, and the output movie file name (e.g., Forced_no_fins_movie).
	This file ties everything together, orchestrating the solver's operation and producing the output movie.


# How It Works

1. Mesh Reading: The Unstructured_Mesh class reads the input mesh file, which contains data about the simulation domain. It stores important information such as node positions, cell connectivity, and boundary conditions, which represent the physical structure (e.g., the heatsink geometry) that will be analyzed.

2. Finite Volume Solver: The Unstructured_FV_Solver class discretizes the governing heat transfer equations (such as the heat conduction equation) using the finite volume method. It calculates the distribution of heat across the mesh. The input heat source, in this case, could be something like the heat emitted from a CPU. The solver integrates these equations over time, iteratively updating the temperature field at each time step based on the mesh geometry and simulation parameters (such as thermal conductivity, heat generation rate, and boundary conditions).

3. Movie Generation: As the simulation progresses, temperature distributions across the mesh are computed at each time step. Images showing these results are generated and compiled into a movie file (e.g., Forced_no_fins_movie.mp4). This movie provides a visual representation of how the temperature changes over time, helping to visualize the heat dissipation process from the CPU or heat source.
