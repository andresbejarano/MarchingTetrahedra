#include "Vector3.h"
#include <math.h>
#include <sstream>

double Vector3::magnitude()
{
	// Calculate the magnitude of the vector and return it
	return sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z));
}

Vector3* Vector3::normalize()
{
	// Get the magnitude of the vector
	double mag = this->magnitude();

	// If it has a magnitude then normalize it
	if (mag > 0)
	{
		return new Vector3(this->x / mag, this->y / mag, this->z / mag);
	}

	// Since the vector has no magnitude then return a zero vector
	return new Vector3(0, 0, 0);
}

Vector3* Vector3::multiply(double s)
{
	// Multiply each component by the scalar
	return new Vector3(this->x * s, this->y * s, this->z * s);
}

Vector3* Vector3::clone()
{
	return new Vector3(this->x, this->y, this->z);
}

Vector3* Vector3::add(Vector3* B)
{
	return new Vector3(this->x + B->x, this->y + B->y, this->z + B->z);
}

Vector3* Vector3::sub(Vector3* B)
{
	return new Vector3(this->x - B->x, this->y - B->y, this->z - B->z);
}

double Vector3::dot(Vector3* B)
{
	return (this->x * B->x) + (this->y * B->y) + (this->z * B->z);
}

Vector3* Vector3::cross(Vector3* B)
{
	// Calculate the cross product values
	double i = (this->y * B->z) - (this->z * B->y);
	double j = (this->z * B->x) - (this->x * B->z);
	double k = (this->x * B->y) - (this->y * B->x);

	// Return the new vector
	return new Vector3(i, j, k);
}

double Vector3::squareDistance(Vector3* B)
{
	return ((B->x - this->x) * (B->x - this->x)) + 
		   ((B->y - this->y) * (B->y - this->y)) + 
		   ((B->z - this->z) * (B->z - this->z));
}

double Vector3::distance(Vector3* B)
{
	return sqrt(squareDistance(B));
}

Vector3* Vector3::interpolate(double t, Vector3* B)
{
	return this->multiply(1.0 - t)->add(B->multiply(t));
}

bool Vector3::equal(Vector3* B)
{
	return (this->x == B->x && this->y == B->y && this->z == B->z);
}

std::string Vector3::toString()
{
	// Convert the vector to a string notation and return it
	std::stringstream ss;
	ss << this->info << ":(" << this->x << ", " << this->y << ", " << this->z << ")";
	return ss.str();
}
