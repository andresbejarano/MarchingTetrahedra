#include "Triangle.h"

bool Triangle::hasVertex(Vector3 * V)
{
	return (v0->equal(V) || v1->equal(V) || v2->equal(V));
}

bool Triangle::isAdjacent(Triangle * T)
{
	// Two triangles are adjacent if they share one common edge, which means they share two vertices
	if (hasVertex(T->v0) && hasVertex(T->v1)) 
	{
		return true;
	}
	else if (hasVertex(T->v1) && hasVertex(T->v2)) 
	{
		return true;
	}
	else if (hasVertex(T->v2) && hasVertex(T->v0)) 
	{
		return true;
	}

	// Since they don't share two vertices then return false
	return false;
}

bool Triangle::isPoint()
{
	return (v0->equal(v1) && v1->equal(v2));
}

Vector3* Triangle::normal()
{
	return v1->sub(v0)->cross(v2->sub(v0));
}

bool Triangle::isCCW(Vector3* V)
{
	// The vertices of the triangle are in CCW if the dot product between the normal of 
	// the triangle and the given vector is positive, which means the angle is positive 
	// and less than 90 degrees
	return (V->dot(normal()) > 0);
}
