#include "mwTPoint3d.hpp"
#include "mwArcFunction.hpp"
#include "CreateSkin.hpp"
//#############################################################################

int main(int argc, char* argv[])
{

	//Point cloud reference point at 0, 0, 0
	cadcam::mwTPoint3d<double> referencePoint( 0., 0., 0. );

	//Number of points in x direction
	const unsigned long nx = 1000;

	//Number of points in y direction
	const unsigned long ny = 100;

	//Number of points in z direction
	const unsigned long nz = 500;

	//Distance between points in the point grid (same fo x, y and z directions)
	const double delta = 1.;

	//Discrete step for move function calculation
	const double deltaT = 0.1;

	//Radius of the sphere
	const double sphereRad = 10.;

	//Accuracity of move (only odd numbers!) bigger=faster
	const double accuracity = 3;

	//Name of the file to write the skin data to
	std::string skinFileName( "test.asc" );

	//Function object to be evaluated
	mwArcFunction func( 0., 1., 150. );

	//Evaluation here
	CreateSkin( referencePoint, nx, ny, nz, sphereRad, func, deltaT,accuracity, delta, skinFileName );

	return 0;
}
