#include<iostream>
#include<fstream>
#include<vector>

using namespace std;

// The following function prints the lattice to file "output.ppm"
void Print_lattice (const vector<int> vlat, const int &vlx, const int &vly, const int &vwidth, const int &vheight, const char* vfilename="output.ppm")
{
  int  i, j, k, l;
  int vw= vwidth/vlx, vh=vheight/vly;
  int r[8], g[8], b[8];

  r[0]= 255; g[0]= 255; b[0]= 255; //white  use 0 in your lattice if you want to color it white
  r[1]= 255; g[1]=   0; b[1]=   0; //red 
  r[2]= 255; g[2]= 128; b[2]=   0; //orange 
  r[3]= 255; g[3]= 255; b[3]=   0; //yellow
  r[4]=   0; g[4]= 255; b[4]=   0; //green
  r[5]=   0; g[5]= 255; b[5]= 255; //magenta
  r[6]=   0; g[6]= 128; b[6]= 255; //light blue
  r[7]=   0; g[7]=   0; b[7]= 255; //blue
  ofstream out (vfilename);

  out << "P3" << endl;
  out << vw*vlx << " " << vh*vly << endl;
  out << "255" << endl;

  for (i=vly-1; i>=0; i--)
    for (j=0; j<vh; j++)
      for (k=0; k<vlx; k++)
      {
        for (l=0; l<vw; l++)
        { out << r[vlat[k+i*vlx]] << " " << g[vlat[k+i*vlx]] << " " << b[vlat[k+i*vlx]] << " ";
        }
      } 
      out << endl;

  out.close ();
}
