
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

	read_mesh: Reads the mesh data from a file (e.g., .msh file format).
	get_node_count: Returns the number of nodes in the mesh.
	get_cell_count: Returns the number of cells in the mesh.
	get_cell_neighbors: Finds the neighboring cells of a given cell for the simulation.
	This class ensures that the mesh data is properly structured and accessible for the solver to work on.

2. Unstructured_FV_Solver
The Unstructured_FV_Solver class is where the core functionality of the finite volume solver resides. It is responsible for setting up and running the simulation, updating the solution at each time step, and generating the output movie.

Key Properties:

	Mesh: An instance of the Unstructured_Mesh class. This provides the solver with the mesh data on which it operates.
	Time Step (dt): The time step used for numerical integration in the simulation.
	Simulation Parameters: Various parameters such as the total simulation time, number of frames to generate, and output file names are defined here.

Key Functions:

	make_movie: This function controls the creation of the output movie. It takes parameters such as the time step, frame interval, and the name of the output movie file. It runs the simulation, updates the solution, and generates images for each frame, which are then compiled into a movie.
	run_simulation: Runs the simulation, solving the physical problem over the mesh. It iterates through time steps, solving for physical quantities (like velocity, pressure, etc.) at each node and cell of the mesh.
	compute_fluxes: Calculates the fluxes across the faces of each cell for the finite volume method.
	The solver uses the unstructured mesh to discretize the governing equations of the physical process and integrates them over time. The results are then stored and visualized.

3. Main Program (main.cpp)
The main program serves as the entry point for the project. It initializes the solver with a mesh file, calls the solver to run the simulation, and generates an animation.

Key Steps:

	Create Solver Instance: An instance of the Unstructured_FV_Solver class is created, with the mesh file (nofins.msh) passed to the constructor.
	Call make_movie: The make_movie function is called to run the simulation and generate the movie. It takes parameters like time step, frame interval, and the output movie file name (e.g., Forced_no_fins_movie).
	This file ties everything together, orchestrating the solver's operation and producing the output movie.


# How It Works

Mesh Reading: The Unstructured_Mesh class reads the input mesh file and stores the mesh data, such as node positions, cell connectivity, and boundary information.
Finite Volume Solver: The Unstructured_FV_Solver class uses the mesh to discretize the governing equations using the finite volume method. The solver integrates over time, updating the solution at each time step based on the mesh geometry and simulation parameters.
Movie Generation: As the simulation progresses, images of the solution at each time step are generated and compiled into a movie file. This provides a visualization of the simulation results.
Example Usage

Once the program is compiled and the necessary files are in place, running the main.cpp will execute the following steps:

Load the mesh file (nofins.msh).
Run the finite volume solver using the mesh and the defined simulation parameters.
Generate a movie (Forced_no_fins_movie.mp4) based on the simulation results.
