#include <iostream>
#include <fstream>


// The following function prints the lattice to file "output.ppm"
void Print_lattice (int *vlat, const int &vlx, const int &vly, const int &vwidth, const int &vheight, const char* vfilename="output.ppm")
{
  const int vw= vwidth / vlx;
  const int vh= vheight / vly;
  int r[5], g[5], b[5];

  r[0]= 255; g[0]= 255; b[0]= 255; //white  use 0 in your lattice if you want to color it white
  r[1]=   0; g[1]= 255; b[1]=   0; //green  use 1 in your lattice if you want to color it green
  r[2]= 255; g[2]=   0; b[2]=   0; //red
  r[3]=   0; g[3]=   0; b[3]=   0; //black
  r[4]=   0; g[4]=   0; b[4]= 255; //blue

  std::ofstream out (vfilename);

  out << "P3" << std::endl;
  out << vw * vlx << " " << vh * vly << std::endl;
  out << "255" << std::endl;

  for (int i=vly-1; i>=0; i--)
    for (int j=0; j<vh; j++)
      for (int k=0; k<vlx; k++)
      {
            for (int l=0; l<vw; l++)
            {
                out << r[vlat[k+i*vlx]] << " " << g[vlat[k+i*vlx]] << " " << b[vlat[k+i*vlx]] << " ";
            }
      } 
      out << std::endl;
  out.close ();
}
