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




bool PointInSphere(cadcam::mwTPoint3d<double> sphereCentr,double sphereRad, cadcam::mwTPoint3d<double> dot) 
{
    return std::pow(sphereCentr.x() - dot.x(), 2) +
        std::pow(sphereCentr.y() - dot.y(), 2) +
        std::pow(sphereCentr.z() - dot.z(), 2) <= std::pow(sphereRad, 2);
}

void CutSphere(mwArcFunction& func, cadcam::mwTPoint3d<double>*** mass,
            const unsigned long nx, const unsigned long ny,
            const unsigned long nz, const double delta, 
            const double deltaT, const double sphereRad,const double accuracity)
{
    //create sphereCenter and path
    cadcam::mwTPoint3d<double> sphereCentr,sphereCentr2,vect,pos,point1,point2;
    std::vector<cadcam::mwTPoint3d<double>> sphereMove;

    double x1, y1, z1, x2, y2, z2;
    
    //loop for function
    for (double dt=func.GetBeginParameter(); dt< func.GetEndParameter(); dt+=deltaT)
    {
        //get center from func
        sphereCentr = func.Evaluate(dt);
        if (dt + deltaT <= func.GetEndParameter())
            sphereCentr2 = func.Evaluate(dt + deltaT);
        else
            sphereCentr2 = func.Evaluate(func.GetEndParameter());
        

        //find vector
        vect = sphereCentr2 - sphereCentr;
        !vect;

        //move sphere from Centr to Centr2 and save the path
        for (pos=sphereCentr; !PointInSphere(sphereCentr2,sphereRad/2,pos); pos += vect*accuracity)
            sphereMove.push_back(pos);

        //define borders for loop
        point1.min(sphereCentr, sphereCentr2);
        point2.max(sphereCentr, sphereCentr2);

        x1 = (point1.x() - sphereRad) / delta;
        y1 = (point1.y() - sphereRad) / delta;
        z1 = (point1.z() - sphereRad) / delta;
        x2 = (point2.x() + sphereRad) / delta;
        y2 = (point2.y() + sphereRad) / delta;
        z2 = (point2.z() + sphereRad) / delta;

        //loop (hide dots)
        for (int i = x1; i < x2; i++) 
        {
            if (i >= nx || i < 0)
                continue;
            for (int k = y1; k < y2; k++) 
            {
                if (k >= ny || k<0)
                    continue;
                for (int j = z1; j < z2 && j < nz && j>0; j++) 
                {
                    if(j >= nz || j < 0)
                        continue;

                    if (PointInSphere(sphereCentr, sphereRad, mass[i][k][j]) ||
                        PointInSphere(sphereCentr2, sphereRad, mass[i][k][j]))
                        mass[i][k][j].setVis(false);

                    for each (cadcam::mwTPoint3d<double> movCentr in sphereMove)
                        if (PointInSphere(movCentr, sphereRad, mass[i][k][j]))
                            mass[i][k][j].setVis(false);

                }
            }
        }
                      
        sphereMove.clear();
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
                mwArcFunction& func, const double deltaT,const double accuracity, //replace mwDiscreteFunction with mwArcFunction! and add accuracity
				const double delta, const std::string &skinFileName )
{
    //create 3d massive
    cadcam::mwTPoint3d<double>*** mass;

    //fill it
    mass=CreateMassive(nx, ny, nz, delta);

    //hide dots
    CutSphere(func,mass, nx, ny, nz, delta, deltaT, sphereRad,accuracity);

    //save result
    SaveSkin(mass, skinFileName, nx, ny, nz);

}
