#ifndef CREATESKIN_HPP
#define CREATESKIN_HPP

#include "mwTPoint3d.hpp"
#include "mwDiscreteFunction.hpp"
#include "mwTPoint3d.hpp"
#include "mwArcFunction.hpp"
#include <fstream>
#include <math.h>
void CreateSkin(const cadcam::mwTPoint3d<double> refPoint,
			const unsigned long nx, const unsigned long ny,
			const unsigned long nz, const double sphereRad,
			mwArcFunction& func, const double deltaT, const double accuracity,
			const double delta, const std::string &skinFileName );

cadcam::mwTPoint3d<double>*** CreateMassive(const unsigned long nx, const unsigned long ny,
	const unsigned long nz, const double delta);

void CutSphere(mwArcFunction& func, cadcam::mwTPoint3d<double>*** mass,
	const unsigned long nx, const unsigned long ny,
	const unsigned long nz, const double delta,
	const double deltaT, const double sphereRad,const double accuracity);

void SaveSkin(cadcam::mwTPoint3d<double>*** mass, std::string skinFileName,
	const unsigned long nx, const unsigned long ny, const unsigned long nz);

#endif /* CREATESKIN_HPP */
