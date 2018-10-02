#include "latticeview.h"
#include <stdlib.h>
#include <stdio.h> // for sprintf()

#define N 10            //Lateral number of cells
#define ImageWidth 1000  //image width
#define ImageHeight 1000 //image height


using namespace std;

int main ()
{
  int lat[N*N];  //create the NxN lattice

  // Several examples: define states and print lattice
  
  /* 
  
  for (int icounter=0; icounter<N*N; icounter++) lat[icounter]=0; //for loop runs through all the cells, setting them to zero
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "allwhite.ppm");

  for (int icounter=0; icounter<N*N; icounter++) lat[icounter]=1; // use 1 for green color
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "allgreen.ppm");

  for (int icounter=0; icounter<N*N; icounter++) lat[icounter]=2; // use 2 for red color
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "allred.ppm");

  for (int icounter=0; icounter<N*N; icounter++) lat[icounter]=3; // use 3 for black color
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "allback.ppm");

  for (int icounter=0; icounter<N*N; icounter++) lat[icounter]=4; // use 4 for blue color
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "allblue.ppm");
  
  for (int icounter=0; icounter<N*N; icounter++)
    lat[icounter]= (int)(drand48()*5.0); // Extracts a random number between 0 and 5.
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "randomcolor.ppm");
  
  */

  // We can try the following example where one row is colored:
  for (int icounter=0; icounter<N; icounter++)
    for (int jcounter=0; jcounter <N; jcounter++)
      lat[icounter+jcounter*N]= jcounter%5;
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "colorrow.ppm");

  // for a series of pictures, you would need VARYING FILENAMES. Use the following lines of code:
  char filename[160];
  int img_no=1;
  sprintf(filename,"test%03d.ppm",img_no); // now same as "colorrow.ppm" - i.e. prints 'colorrow' as 'test001'
  Print_lattice (lat, N, N, ImageWidth, ImageHeight,filename);
  
  // Note, that you should have the correct order of the files. (if it is more than 10, it should be test01.png ....). You may have to use: sprintf(filename_png,"test%03d.png",img_no); 

  for (int icounter=0; icounter<N; icounter++)
    for (int jcounter=0; jcounter <N; jcounter++)
      lat[icounter+jcounter*N]= icounter%5;
  Print_lattice (lat, N, N, ImageWidth, ImageHeight, "colorcolumn.ppm");
  
  // the next image of the series
  img_no++;
  sprintf(filename,"test%03d.ppm",img_no); // now same as "colorcolumn.ppm" - i.e. prints 'colorcolumn' as 'test002'
  Print_lattice (lat, N, N, ImageWidth, ImageHeight,filename);
  
  
  // In order to save memory you can directly convert the ppm's to any other format like, e.g. png. Use the following code to convert ppm to another format (e.g. png, could be jpg, etc.)
  char cmd[160],filename_png[160]; 
  sprintf(filename_png,"test%03d.png",img_no);
  sprintf(cmd,"convert %s %s; rm -f %s",filename,filename_png,filename); // rm automatically deletes the ppm-file to save memory (removes test002.ppm)
  
  // sprintf(cmd,"convert %s %s",filename,filename_png,filename); // only converts
  system(cmd);
  
  // From a sequence of pictures you can easily create a movie, i.e. an animated gif by simply using the "convert" command, i.e. the following command line in your terminal/shell: convert test*.png test.gif
}
