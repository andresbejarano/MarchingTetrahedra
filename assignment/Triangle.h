#pragma once

#include "Vector3.h"

class Triangle 
{
	public:

		// The vertices of the triangle
		Vector3* v0;
		Vector3* v1;
		Vector3* v2;

		/*
			Constructor of the class
		*/
		Triangle(Vector3* v0, Vector3* v1, Vector3* v2) : v0(v0), v1(v1), v2(v2) {}

		/*
			Indicates whether the triangle has a vertex with the same (x, y, z) coordinate values
			@param V A vertex
			@return true if the triangle has a vertex with the same (x, y, z) values, otherwise false
		*/
		bool hasVertex(Vector3* V);

		/*
			Indicates whether the triangle is adjacent (shares an edge) with triangle T
			@param T A triangle
			@return true if both triangles share an edge, otherwise false
		*/
		bool isAdjacent(Triangle* T);

		/*
			Indicates whether the triangle's vertices are the same point
			@return true if the vertices are the same, otherwise false
		*/
		bool isPoint();

		/*
			Returns the normal vector of the triangle using the vertices
			@return a Vector3 pointer to the object representing the normal
		*/
		Vector3* normal();

		/*
			Indicates whether the triangle has counterclockwise direction with respect of a given vector
			@param V A vector
			@return true if the triangle vertices are in counterclockwise direction with respect of the 
			given vector
		*/
		bool isCCW(Vector3* V);
};
