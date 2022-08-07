#include "CreateSkin.hpp"


cadcam::mwTPoint3d<double>*** CreateMassive( const unsigned long nx, const unsigned long ny,
					                        const unsigned long nz, const double delta) 
{
    //add variables for each dot`s dimention and create 3d massive of dots
	float dx = 0, dy = 0, dz = 0;
	cadcam::mwTPoint3d<double>*** mass = new cadcam::mwTPoint3d<double>**[nx];

    //init the massive
	for (int i = 0; i < nx; i++)
	{
		mass[i] = new cadcam::mwTPoint3d<double>*[ny];
		for (int j = 0; j < ny; j++)
			mass[i][j] = new cadcam::mwTPoint3d<double>[nz];
	}

    //fill it with data
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                mass[i][j][k] = { dx, dy, dz };
                dz += delta;
            }
            dy += delta;
            dz = 0;
        }
        dx += delta;
        dy = 0;
    }
    return mass;
}


void CutSphere(mwArcFunction& func, cadcam::mwTPoint3d<double>*** mass,
            const unsigned long nx, const unsigned long ny,
            const unsigned long nz, const double delta, 
            const double deltaT, const double sphereRad)
{
    //create sphereCenter
    cadcam::mwTPoint3d<double> sphereCentr;
    int x1, y1, z1, x2, y2, z2;

    //loop for function
    for (double dt=func.GetBeginParameter(); dt< func.GetEndParameter(); dt+=deltaT)
    {
        //get center from func
        sphereCentr = {func.Evaluate(dt)};

        //define borders for loop
        x1 = sphereCentr.x() - sphereRad / delta;
        y1 = sphereCentr.y() - sphereRad / delta;
        z1 = sphereCentr.z() - sphereRad / delta;
        x2 = sphereCentr.x() + sphereRad / delta;
        y2 = sphereCentr.y() + sphereRad / delta;
        z2 = sphereCentr.z() + sphereRad / delta;

        for (int i = x1; i < x2 && i < nx; i++)
            for (int k = y1; k < y2 && k < ny; k++)
                for (int j = z1; j < z2 && j < nz; j++)
                {
                    if (std::pow(sphereCentr.x() - mass[i][k][j].x(), 2) +
                        std::pow(sphereCentr.y() - mass[i][k][j].y(), 2) +
                        std::pow(sphereCentr.z() - mass[i][k][j].z(), 2) <= std::pow(sphereRad, 2))
                        mass[i][k][j].setVis(false);
                }
    }
    
}




void SaveSkin( cadcam::mwTPoint3d<double>*** mass, std::string skinFileName, 
              const unsigned long nx, const unsigned long ny,const unsigned long nz)
{
    //create data buffer and save stream
    std::string data = "";
    std::ofstream saveFile(skinFileName);

    int yl;

    for (int i = 0; i < nx; i++)
    {

        for (int k = 0; k < nz; k++)
        {
            yl = 0;

            //select the highest dot 
            while (yl < ny - 1 && mass[i][yl + 1][k].vis())
            {
                mass[i][yl++][k].setVis(false);
            }

            //add selected dot to buffer
            data += std::to_string(mass[i][yl][k].x()) + ' ' + std::to_string(mass[i][yl][k].y()) + ' ' + std::to_string(mass[i][yl][k].z()) + '\n';
        }
    }

    //put data in file and save it 
    saveFile << data;
    saveFile.close();
}

 

void CreateSkin( const cadcam::mwTPoint3d<double> refPoint, 
				const unsigned long nx, const unsigned long ny, 
				const unsigned long nz, const double sphereRad, 
                mwArcFunction& func, const double deltaT, //replace mwDiscreteFunction with mwArcFunction!
				const double delta, const std::string &skinFileName )
{
    //create 3d massive
    cadcam::mwTPoint3d<double>*** mass;

    //fill it
    mass=CreateMassive(nx, ny, nz, delta);

    //hide dots
    CutSphere(func,mass, nx, ny, nz, delta, deltaT, sphereRad);

    //save result
    SaveSkin(mass, skinFileName, nx, ny, nz);

}
