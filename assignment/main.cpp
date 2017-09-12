#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle, vtkRenderingFreeType, vtkRenderingOpenGL2)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2)

#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTriangle.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkStructuredPoints.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkIdList.h>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <queue>
#include "Vector3.h"
#include "Triangle.h"


/*	
	+------------------------------+
	| GLOBAL VARIABLES AND OBJECTS |
	+------------------------------+
*/

// Define the objects where the points and triangles will be stored
vtkSmartPointer<vtkPoints> isosurfacePoints = vtkSmartPointer<vtkPoints>::New();
vtkSmartPointer<vtkCellArray> isosurfaceTriangles = vtkSmartPointer<vtkCellArray>::New();

// Keep record of the location for each triangle vertex in the isosurface points object
// NOTE: This is done for avoiding additional/redundant points in the final mesh
std::unordered_map<std::string, vtkIdType>* pointIndex = new std::unordered_map<std::string, vtkIdType>();

// Keep record of the IDs given to the triangles
// NOTE: We need them for finding the shells quickly. The bool variable indicates whether the triangle has already been 
// visited during shell discovery
std::unordered_map<vtkIdType, int>* triangleIndex = new std::unordered_map<vtkIdType, int>();
std::unordered_map<vtkIdType, vtkIdType>* triangleV0Index = new std::unordered_map<vtkIdType, vtkIdType>();
std::unordered_map<vtkIdType, vtkIdType>* triangleV1Index = new std::unordered_map<vtkIdType, vtkIdType>();
std::unordered_map<vtkIdType, vtkIdType>* triangleV2Index = new std::unordered_map<vtkIdType, vtkIdType>();

vtkSmartPointer<vtkPolyData> isosurfacePolydata;



/*
	+-----------+
	| FUNCTIONS |
	+-----------+
*/


/*
	Returns the parameter t such that it gives the linear interpolation value f between f1 and f2
*/
double inverseLinearInterpolation(double f, double f1, double f2) 
{
	return (f - f1) / (f2 - f1);
}


/*
	Builds the triangle corresponding to the indicated vertices and given isovalue
	@param V0 Vertex that has a different value from the other three vertices. If this is positive 
	the other vertices must be negative, and viceversa
	@param V1 The next vertex in counter clockwise order
	@param V2 The next next vertex in counter clockwise order
	@param V3 The vertex that is not on the plane where V0, V1 and V2 are
	@param isovalue The isovalue to be considered for inverse interpolation
	@return The pointer to the generated triangle
*/
Triangle* buildTriangle(Vector3* V0, Vector3* V1, Vector3* V2, Vector3* V3, double isovalue)
{
	// Get the t parameters for the intersected values in the edges
	double t01 = inverseLinearInterpolation(isovalue, V0->info, V1->info);
	double t02 = inverseLinearInterpolation(isovalue, V0->info, V2->info);
	double t03 = inverseLinearInterpolation(isovalue, V0->info, V3->info);

	// Interpolate values and get the triangle vertices
	Vector3* T0 = V0->interpolate(t01, V1);
	Vector3* T1 = V0->interpolate(t03, V3);
	Vector3* T2 = V0->interpolate(t02, V2);

	// Generate the triangle
	Triangle* T = new Triangle(T0, T1, T2);

	// If it happens to be a degenerate triangle (all the vertices are in the same point) then 
	// return NULL, Otherwise return the triangle
	return (T->isPoint()) ? NULL : T;
}

/*
	Builds the two triangles corresponding to the indicated vertices and given isovalue
	@param V0 One of the vertices that has a different value from the other two vertices. If this is positive
	the other vertices must be negative, and viceversa
	@param V1 The other vertex with the same sign value as V0
	@param V2 The next next vertex in counter clockwise order
	@param V3 The vertex that is not on the plane where V0, V1 and V2 are
	@param isovalue The isovalue to be considered for inverse interpolation
	@return The pointer to a vector with the generated triangles
*/
std::vector<Triangle*>* buildTriangles(Vector3* V0, Vector3* V1, Vector3* V2, Vector3* V3, double isovalue)
{
	// Get the t parameters for the intersected values in the edges
	double t02 = inverseLinearInterpolation(isovalue, V0->info, V2->info);
	double t03 = inverseLinearInterpolation(isovalue, V0->info, V3->info);
	double t12 = inverseLinearInterpolation(isovalue, V1->info, V2->info);
	double t13 = inverseLinearInterpolation(isovalue, V1->info, V3->info);

	Vector3* T0 = V0->interpolate(t02, V2);
	Vector3* T1 = V1->interpolate(t12, V2);
	Vector3* T2 = V1->interpolate(t13, V3);
	Vector3* T3 = V0->interpolate(t03, V3);

	// Generate the triangles and store them in a vector
	Triangle* triangle1 = new Triangle(T0, T1, T2);
	Triangle* triangle2 = new Triangle(T2, T3, T0);
	std::vector<Triangle*>* triangles = new std::vector<Triangle*>();
	
	// If the triangles are not degenerate (their vertices are the same point) then add them to the triangles vector
	if (!triangle1->isPoint())
	{
		triangles->push_back(triangle1);
	}

	if (!triangle2->isPoint())
	{
		triangles->push_back(triangle2);
	}
	
	// Returns the triangles
	return triangles;
}


/*
	March the tetrahedron given by its vertices and returns the generated triangles
	@param V0
	@param V1
	@param V2
	@param V3
	@param isovalue
	@return A pointer to a vector with the generated triangles
*/
std::vector<Triangle*>* marchTetrahedra(Vector3* V0, Vector3* V1, Vector3* V2, Vector3* V3, double isovalue)
{
	// Initialize the vector where the triangles will be stored
	std::vector<Triangle*>* triangles = new std::vector<Triangle*>();

	// Define the signs of each vertex of the tetrahedron
	int s0 = (V0->info > isovalue) ? 1 : 0;
	int s1 = (V1->info > isovalue) ? 1 : 0;
	int s2 = (V2->info > isovalue) ? 1 : 0;
	int s3 = (V3->info > isovalue) ? 1 : 0;

	// Process each case
	if (s0 == 0 && s1 == 0 && s2 == 0 && s3 == 0) 
	{
		// CASE 0: All negative => No triangles
	}
	else if (s0 == 1 && s1 == 0 && s2 == 0 && s3 == 0)
	{
		// Case 1: 1000
		// One positive vertex, all others negative => One triangle
		// Vertices are in edges 0-1, 0-2 and 0-3

		// Generate the triangle and add it only if it is not NULL
		// NOTE: A triangle is NULL if it is a degenerated triangle
		Triangle* triangle = buildTriangle(V0, V1, V2, V3, isovalue);
		if (triangle != NULL) 
		{
			triangles->push_back(triangle);
		}
	}
	else if (s0 == 0 && s1 == 1 && s2 == 0 && s3 == 0)
	{
		// Case 2: 0100
		// One positive vertex, all others negative => One triangle
		// Vertices are in edges 1-2, 1-3 and 1-0

		// Generate the triangle and add it only if it is not NULL
		// NOTE: A triangle is NULL if it is a degenerated triangle
		Triangle* triangle = buildTriangle(V1, V2, V0, V3, isovalue);
		if (triangle != NULL)
		{
			triangles->push_back(triangle);
		}
	}
	else if (s0 == 1 && s1 == 1 && s2 == 0 && s3 == 0)
	{
		std::vector<Triangle*>* t = buildTriangles(V0, V1, V2, V3, isovalue);
		for (auto it = t->begin(); it != t->end(); ++it) 
		{
			triangles->push_back(*it);
		}
	}
	else if (s0 == 0 && s1 == 0 && s2 == 1 && s3 == 0)
	{
		// Case 4: 0010
		// One positive vertex, all others negative => One triangle
		// Vertices are in edges 2-0, 2-1 and 2-3

		// Generate the triangle and add it only if it is not NULL
		// NOTE: A triangle is NULL if it is a degenerated triangle
		Triangle* triangle = buildTriangle(V2, V0, V1, V3, isovalue);
		if (triangle != NULL)
		{
			triangles->push_back(triangle);
		}
	}
	else if (s0 == 1 && s1 == 0 && s2 == 1 && s3 == 0)
	{
		std::vector<Triangle*>* t = buildTriangles(V0, V2, V1, V3, isovalue);
		for (auto it = t->begin(); it != t->end(); ++it)
		{
			triangles->push_back(*it);
		}
	}
	else if (s0 == 0 && s1 == 1 && s2 == 1 && s3 == 0)
	{
		std::vector<Triangle*>* t = buildTriangles(V1, V2, V0, V3, isovalue);
		for (auto it = t->begin(); it != t->end(); ++it)
		{
			triangles->push_back(*it);
		}
		
	}
	else if (s0 == 1 && s1 == 1 && s2 == 1 && s3 == 0)
	{
		// Case 7: 1110
		// One negative vertex, all others positive => One triangle
		// Vertices are in edges 3-2, 3-0 and 3-1

		// Generate the triangle and add it only if it is not NULL
		// NOTE: A triangle is NULL if it is a degenerated triangle
		Triangle* triangle = buildTriangle(V3, V0, V2, V1, isovalue);
		if (triangle != NULL)
		{
			triangles->push_back(triangle);
		}
	}
	else if (s0 == 0 && s1 == 0 && s2 == 0 && s3 == 1)
	{
		// Case 8: 0001
		// One positive vertex, all others negative => One triangle
		// Vertices are in edges 3-2, 3-0 and 3-1

		// Generate the triangle and add it only if it is not NULL
		// NOTE: A triangle is NULL if it is a degenerated triangle
		Triangle* triangle = buildTriangle(V3, V2, V1, V0, isovalue);
		if (triangle != NULL)
		{
			triangles->push_back(triangle);
		}
	}
	else if (s0 == 1 && s1 == 0 && s2 == 0 && s3 == 1)
	{
		std::vector<Triangle*>* t = buildTriangles(V0, V3, V1, V2, isovalue);
		for (auto it = t->begin(); it != t->end(); ++it)
		{
			triangles->push_back(*it);
		}
	}
	else if (s0 == 0 && s1 == 1 && s2 == 0 && s3 == 1)
	{
		std::vector<Triangle*>* t = buildTriangles(V1, V3, V0, V2, isovalue);
		for (auto it = t->begin(); it != t->end(); ++it)
		{
			triangles->push_back(*it);
		}
	}
	else if (s0 == 1 && s1 == 1 && s2 == 0 && s3 == 1)
	{
		// Case 11: 1101
		// One negative vertex, all others positive => One triangle
		// Vertices are in edges 2-0, 2-3 and 2-1

		// Generate the triangle and add it only if it is not NULL
		// NOTE: A triangle is NULL if it is a degenerated triangle
		Triangle* triangle = buildTriangle(V2, V0, V1, V3, isovalue);
		if (triangle != NULL)
		{
			triangles->push_back(triangle);
		}
	}
	else if (s0 == 0 && s1 == 0 && s2 == 1 && s3 == 1)
	{
		std::vector<Triangle*>* t = buildTriangles(V2, V3, V0, V1, isovalue);
		for (auto it = t->begin(); it != t->end(); ++it)
		{
			triangles->push_back(*it);
		}
	}
	else if (s0 == 1 && s1 == 0 && s2 == 1 && s3 == 1)
	{
		// Case 13: 1011
		// One negative vertex, all others positive => One triangle
		// Vertices are in edges 1-2, 1-3 and 1-0

		// Generate the triangle and add it only if it is not NULL
		// NOTE: A triangle is NULL if it is a degenerated triangle
		Triangle* triangle = buildTriangle(V1, V2, V0, V3, isovalue);
		if (triangle != NULL)
		{
			triangles->push_back(triangle);
		}
		
	}
	else if (s0 == 0 && s1 == 1 && s2 == 1 && s3 == 1)
	{
		// Case 14: 0111
		// One negative vertex, all others positive => One triangle
		// Vertices are in edges 0-1, 0-3 and 0-2
		
		// Generate the triangle and add it only if it is not NULL
		// NOTE: A triangle is NULL if it is a degenerated triangle
		Triangle* triangle = buildTriangle(V0, V1, V2, V3, isovalue);
		if (triangle != NULL)
		{
			triangles->push_back(triangle);
		}
	}
	else if (s0 == 1 && s1 == 1 && s2 == 1 && s3 == 1)
	{
		// CASE 15: All positive => No triangles
	}
	else 
	{
		// Something weird is occuring here!
	}

	// Return the generated triangles
	return triangles;
}

/*
	Run the marching tetrahedra using the information of the cell vertices and the isovalue
	@param v0
	@param v1
	@param v2
	@param v3
	@param v4
	@param v5
	@param v6
	@param isovalue
	@return The pointer to the vector containing the triangles from the cell
*/
std::vector<Triangle*>* marchCellTetrahedra(Vector3* v0, Vector3* v1, Vector3* v2, Vector3* v3, Vector3* v4, Vector3* v5, Vector3* v6, Vector3* v7, double isovalue) 
{
	// Initialize the vector where the triangles will be stored
	std::vector<Triangle*>* triangles = new std::vector<Triangle*>();

	// Run the march algorithm on each tetrahedra and keep the generated triangles
	std::vector<Triangle*>* tetrahedra1Triangles = marchTetrahedra(v0, v1, v3, v5, isovalue);
	std::vector<Triangle*>* tetrahedra2Triangles = marchTetrahedra(v1, v2, v3, v5, isovalue);
	std::vector<Triangle*>* tetrahedra3Triangles = marchTetrahedra(v0, v3, v4, v5, isovalue);
	std::vector<Triangle*>* tetrahedra4Triangles = marchTetrahedra(v2, v3, v5, v6, isovalue);
	std::vector<Triangle*>* tetrahedra5Triangles = marchTetrahedra(v3, v4, v5, v7, isovalue);
	std::vector<Triangle*>* tetrahedra6Triangles = marchTetrahedra(v3, v5, v6, v7, isovalue);

	// Merge the tetrahedra vectors into the general one
	for (auto it = tetrahedra1Triangles->begin(); it != tetrahedra1Triangles->end(); ++it) 
	{
		triangles->push_back(*it);
	}

	for (auto it = tetrahedra2Triangles->begin(); it != tetrahedra2Triangles->end(); ++it)
	{
		triangles->push_back(*it);
	}

	for (auto it = tetrahedra3Triangles->begin(); it != tetrahedra3Triangles->end(); ++it)
	{
		triangles->push_back(*it);
	}

	for (auto it = tetrahedra4Triangles->begin(); it != tetrahedra4Triangles->end(); ++it)
	{
		triangles->push_back(*it);
	}

	for (auto it = tetrahedra5Triangles->begin(); it != tetrahedra5Triangles->end(); ++it)
	{
		triangles->push_back(*it);
	}

	for (auto it = tetrahedra6Triangles->begin(); it != tetrahedra6Triangles->end(); ++it)
	{
		triangles->push_back(*it);
	}

	// Return the triangles
	return triangles;
}


/*
	Traverses from the given triangle Id all the neighbors and adds them to a set. Recursively traverses 
	the mesh until all possible triangles are visited
	@param triangleId The id of the triangle to be considered
	@param sheelTriangles The set of triangles than have been found to belong to the same shell
	@return
*/
std::set<vtkIdType>* findShell(vtkIdType triangleId)
{
	// Initialize the set where the shell triangles will be stored
	std::set<vtkIdType>* shellTriangles = new std::set<vtkIdType>();

	// Initialize the queue and insert the given triangle id
	std::queue<vtkIdType>* queue = new std::queue<vtkIdType>();
	queue->push(triangleId);

	// Repeat while the queue is not empty
	while (!queue->empty()) 
	{
		// Remove the first element from the queue
		vtkIdType id = queue->front();
		//std::cout << "From queue: " << id << std::endl;
		queue->pop();

		// Insert the triangle id into the shell triangles
		shellTriangles->insert(id);
		//std::cout << "Shell triangles: " << shellTriangles->size() << std::endl;

		// Indicate the current triangle is visited
		triangleIndex->at(id) = 1;
		//std::cout << "Should be true but it is " << triangleIndex->at(id) << std::endl;

		// Get the id of the points that compose the current triangle
		vtkSmartPointer<vtkIdList> trianglePointsId = vtkSmartPointer<vtkIdList>::New();
		isosurfacePolydata->GetCellPoints(id, trianglePointsId);

		// Traverse through the points of the triangles
		for (vtkIdType i = 0; i < trianglePointsId->GetNumberOfIds(); i++) 
		{
			// Define a list containing the current triangle point
			vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
			idList->InsertNextId(trianglePointsId->GetId(i));

			// Initialize the list where the neighbor triangles will be stored
			vtkSmartPointer<vtkIdList> neighborTriangles = vtkSmartPointer<vtkIdList>::New();

			// Get the id of the neighbors triangles
			isosurfacePolydata->GetCellNeighbors(id, idList, neighborTriangles);

			// Traverse through the neighbor triangles
			for (vtkIdType j = 0; j < neighborTriangles->GetNumberOfIds(); j++) 
			{
				// Get the id of the current neighbor triangle
				vtkIdType neighborId = neighborTriangles->GetId(j);

				// If the current neighbor triangle is not visited then add it to the queue
				if (triangleIndex->find(neighborId) == triangleIndex->end())
				{
					//std::cout << neighborId << "Doesn't exist" << std::endl;
				}
				else if(triangleIndex->at(neighborId) == 1)
				{
					//std::cout << neighborId << " already visited" << std::endl;
				}
				else 
				{
					//std::cout << neighborId << " to the queue" << std::endl;
					triangleIndex->at(neighborId) = 1;
					queue->push(neighborId);
				}
			}
		}

		//std::cout << queue->size() << std::endl;
	}

	// Return the set with the shell triangles
	return shellTriangles;
	
}


/*
	The main function
*/
int main(int argc, char** argv)
{
	// Read the isovalue
	double isovalue = 0.0;
	std::cout << "Isovalue = ";
	std::cin >> isovalue;
	//double isovalue = 140;

	// Define the structuerd point reader and read the file
	vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
	reader->SetFileName("heart.vtk");
	reader->Update();

	// Get the points, the dimensions and the spacing of the grid
	// NOTE: We need the spacing since the original vtk file is not a squared grid but a rectangular grid
	vtkSmartPointer<vtkStructuredPoints> structuredPoints = reader->GetOutput();
	int* dimensions = structuredPoints->GetDimensions();
	double* spacing = structuredPoints->GetSpacing();

	// Generate the 3D data structure where the triangles of each cell will be stored
	// NOTE: This is done for finding the shells
	std::vector<Triangle*>**** triangles = new std::vector<Triangle*>***[dimensions[0] - 1];
	for (int i = 0; i < dimensions[0] - 1; i += 1) 
	{
		triangles[i] = new std::vector<Triangle*>**[dimensions[1] - 1];

		for (int j = 0; j < dimensions[1] - 1; j += 1) 
		{
			triangles[i][j] = new std::vector<Triangle*>*[dimensions[2] - 1];
		}
	}

	// Traverse the points
	// NOTE: Actually, we are traversing by cells. But since the information is given in points then let's keep such approach
	for (int z = 0; z < dimensions[2] - 1; z += 1) 
	{
		for (int y = 0; y < dimensions[1] - 1; y += 1) 
		{
			for (int x = 0; x < dimensions[0] - 1; x += 1) 
			{
				// Get the values of the vertices of the current cell and apply the spacing to the vertices
				Vector3* v0 = new Vector3(x, y, z);
				v0->info = static_cast<float*>(structuredPoints->GetScalarPointer(v0->x, v0->y, v0->z))[0];
				v0->x *= spacing[0];
				v0->y *= spacing[1];
				v0->z *= spacing[2];

				Vector3* v1 = new Vector3(x + 1, y, z);
				v1->info = static_cast<float*>(structuredPoints->GetScalarPointer(v1->x, v1->y, v1->z))[0];
				v1->x *= spacing[0];
				v1->y *= spacing[1];
				v1->z *= spacing[2];

				Vector3* v2 = new Vector3(x + 1, y, z + 1);
				v2->info = static_cast<float*>(structuredPoints->GetScalarPointer(v2->x, v2->y, v2->z))[0];
				v2->x *= spacing[0];
				v2->y *= spacing[1];
				v2->z *= spacing[2];

				Vector3* v3 = new Vector3(x, y, z + 1);
				v3->info = static_cast<float*>(structuredPoints->GetScalarPointer(v3->x, v3->y, v3->z))[0];
				v3->x *= spacing[0];
				v3->y *= spacing[1];
				v3->z *= spacing[2];

				Vector3* v4 = new Vector3(x, y + 1, z);
				v4->info = static_cast<float*>(structuredPoints->GetScalarPointer(v4->x, v4->y, v4->z))[0];
				v4->x *= spacing[0];
				v4->y *= spacing[1];
				v4->z *= spacing[2];

				Vector3* v5 = new Vector3(x + 1, y + 1, z);
				v5->info = static_cast<float*>(structuredPoints->GetScalarPointer(v5->x, v5->y, v5->z))[0];
				v5->x *= spacing[0];
				v5->y *= spacing[1];
				v5->z *= spacing[2];

				Vector3* v6 = new Vector3(x + 1, y + 1, z + 1);
				v6->info = static_cast<float*>(structuredPoints->GetScalarPointer(v6->x, v6->y, v6->z))[0];
				v6->x *= spacing[0];
				v6->y *= spacing[1];
				v6->z *= spacing[2];

				Vector3* v7 = new Vector3(x, y + 1, z + 1);
				v7->info = static_cast<float*>(structuredPoints->GetScalarPointer(v7->x, v7->y, v7->z))[0];
				v7->x *= spacing[0];
				v7->y *= spacing[1];
				v7->z *= spacing[2];

				// March the cell's tetrahedra and store the triangles in the respective cell
				triangles[x][y][z] = marchCellTetrahedra(v0, v1, v2, v3, v4, v5, v6, v7, isovalue);

				// Traverse the obtained triangles and add them to the isosurface
				// NOTE: Before adding each triangle it is checked whether any vertex is already in the isosurface
				for (auto it = triangles[x][y][z]->begin(); it != triangles[x][y][z]->end(); ++it) 
				{
					// Get the pointer to the current triangle
					Triangle* triangle = (*it);

					// Initialize the location of each triangle vertex in the isosurface points object
					vtkIdType v0Position = -1;
					vtkIdType v1Position = -1;
					vtkIdType v2Position = -1;

					// If the vertex v0 is not found in the point index structure then add it to the isosurface 
					// points and keep its index for later triangle generation
					std::unordered_map<std::string, vtkIdType>::const_iterator got = pointIndex->find(triangle->v0->toString());
					if (got == pointIndex->end()) 
					{
						// Add the point to the isosurface points object
						v0Position = isosurfacePoints->InsertNextPoint(triangle->v0->x, triangle->v0->y, triangle->v0->z);

						// Store the location of the current vertex in the isosurface points object
						pointIndex->insert({ triangle->v0->toString() , v0Position });
					}
					else 
					{
						// Get the location of the vertex in the isosurface points object
						v0Position = got->second;
					}

					// If the vertex v1 is not found in the point index structure then add it to the isosurface 
					// points and keep its index for later triangle generation
					got = pointIndex->find(triangle->v1->toString());
					if (got == pointIndex->end())
					{
						// Add the point to the isosurface points object
						v1Position = isosurfacePoints->InsertNextPoint(triangle->v1->x, triangle->v1->y, triangle->v1->z);

						// Store the location of the current vertex in the isosurface points object
						pointIndex->insert({ triangle->v1->toString() , v1Position });
					}
					else
					{
						// Get the location of the vertex in the isosurface points object
						v1Position = got->second;
					}

					// If the vertex v2 is not found in the point index structure then add it to the isosurface 
					// points and keep its index for later triangle generation
					got = pointIndex->find(triangle->v2->toString());
					if (got == pointIndex->end())
					{
						// Add the point to the isosurface points object
						v2Position = isosurfacePoints->InsertNextPoint(triangle->v2->x, triangle->v2->y, triangle->v2->z);

						// Store the location of the current vertex in the isosurface points object
						pointIndex->insert({ triangle->v2->toString() , v2Position });
					}
					else
					{
						// Get the location of the vertex in the isosurface points object
						v2Position = got->second;
					}

					// Generate the triangle and add it to the triangles object. Save the triangle id as well
					vtkSmartPointer<vtkTriangle> isosurfaceTriangle = vtkSmartPointer<vtkTriangle>::New();
					isosurfaceTriangle->GetPointIds()->SetId(0, v0Position);
					isosurfaceTriangle->GetPointIds()->SetId(1, v1Position);
					isosurfaceTriangle->GetPointIds()->SetId(2, v2Position);
					vtkIdType triangleId = isosurfaceTriangles->InsertNextCell(isosurfaceTriangle);
					triangleIndex->insert({ triangleId, 0 });
					triangleV0Index->insert({ triangleId, v0Position });
					triangleV1Index->insert({ triangleId, v1Position });
					triangleV2Index->insert({ triangleId, v2Position });
				}
			}
		}
	}

	// Define the polydata object for the isosurface
	isosurfacePolydata = vtkSmartPointer<vtkPolyData>::New();
	isosurfacePolydata->SetPoints(isosurfacePoints);
	isosurfacePolydata->SetPolys(isosurfaceTriangles);

	// Define the mapper for the isosurface
	vtkSmartPointer<vtkPolyDataMapper> mapper;

	// Ask the user for a shell finding option
	int option;
	std::cout << "Type 1 for vtkPolyDataConnectivityFilter or 2 for manual shell traversal: ";
	std::cin >> option;

	if (option == 1) 
	{
		/*
		+----------------------------+
		| OPTION 1:                  |
		| Let vtk to find the shells |
		+----------------------------+
		*/

		// Define a filter for finding the shells
		vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
		confilter->SetInputData(isosurfacePolydata);
		confilter->SetExtractionModeToAllRegions();
		confilter->ColorRegionsOn();

		// Define the mapper for the polydata
		mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(confilter->GetOutputPort());
	}
	else 
	{
		/*
		+--------------------------+
		| OPTION 2:                |
		| Find the shells manually |
		+--------------------------+
		*/

		bool allVisited = false;
		int numShells = 0;
		std::vector<std::set<vtkIdType>*>* foundMeshes = new std::vector<std::set<vtkIdType>*>();

		// While all triangles haven't been visited
		while (!allVisited)
		{
			// Assume all triangles have been visited
			allVisited = true;

			// Traverse through the triangles
			for (auto it = triangleIndex->begin(); it != triangleIndex->end(); ++it)
			{
				// If the current triangle hasn't been visited
				if ((*it).second == 0)
				{
					// Indicate a triangle was found that hasn't been visited
					allVisited = false;

					// Find the shell starting on the current triangle
					vtkIdType triangleId = (*it).first;
					std::set<vtkIdType>* shellTriangles = findShell(triangleId);
					foundMeshes->push_back(shellTriangles);

					std::cout << "Shell found! Number of triangles: " << shellTriangles->size() << std::endl;

					// Increment the number of found meshes by one
					numShells += 1;

					break;
				}
			}
		}


		int selectedShell = 0;

		// If more than one shell was found then ask which one is to be visualized
		if (numShells > 1)
		{
			std::cout << "It was found " << numShells << " shells in the isosurface." << std::endl;
			std::cout << "Indicate which shell do you want to visualize (0 for the complete isosurface): ";
			std::cin >> selectedShell;
		}

		// If the number of shells found is one then show it directly
		if (selectedShell == 0)
		{
			// Load the shell (it is the same as the complete mesh
			mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(isosurfacePolydata);
		}
		else
		{
			// Reduce the index of the selected shell
			selectedShell -= 1;

			// Get the triangle ids from the selected shell
			std::set<vtkIdType>* shellTriangleIds = foundMeshes->at(selectedShell);
			vtkSmartPointer<vtkCellArray> shellTriangles = vtkSmartPointer<vtkCellArray>::New();

			for (auto it = shellTriangleIds->begin(); it != shellTriangleIds->end(); ++it)
			{
				// Get the id of the triangle
				vtkIdType id = (*it);

				// Generate the triangle and add it to the triangles object. Save the triangle id as well
				vtkSmartPointer<vtkTriangle> isosurfaceTriangle = vtkSmartPointer<vtkTriangle>::New();
				isosurfaceTriangle->GetPointIds()->SetId(0, triangleV0Index->at(id));
				isosurfaceTriangle->GetPointIds()->SetId(1, triangleV1Index->at(id));
				isosurfaceTriangle->GetPointIds()->SetId(2, triangleV2Index->at(id));
				shellTriangles->InsertNextCell(isosurfaceTriangle);
			}

			// Define the polydata for the shell
			isosurfacePolydata = vtkSmartPointer<vtkPolyData>::New();
			isosurfacePolydata->SetPoints(isosurfacePoints);
			isosurfacePolydata->SetPolys(shellTriangles);

			// Load the shell into the mapper
			mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(isosurfacePolydata);
		}
	}

	// Define the actor for the mapper
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	
	// Define the renderer and add the actor for the isosurface
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->SetBackground(0.1, 0.1, 0.1);

	// Define the renderer window
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindow->SetSize(800, 800);
	renderWindow->Render();

	// Define the renderer window interactor and start the program
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}
