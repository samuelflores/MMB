/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

/* fread example: read a complete file */
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "mrc.h"

using namespace std;

int main () {
  FILE * pFile;
  long lSize;
  char * buffer;
  size_t result;
  // was:
  //pFile = fopen ( "myfile.bin" , "rb" );
  //pFile = fopen ( "/usr/local/ccp4-6.1.2/examples/data/1rxf.mtz" , "rb"  );
  pFile = fopen ( "/usr/local/ccp4-6.1.2/examples/data/SSADinsulin.mtz" , "rb"  );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
  fseek (pFile , 0 , SEEK_SET);
  
  
  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile); // size in bytes
  rewind (pFile);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*1024 ); //sizeof(char) = 1
  mrcH header;

  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  //char  temp[4];
  //fread (&temp,1,4,pFile); 
  result = fread (&header,1,1024 ,pFile); 
  //You need to deal with the optional header extension:
  char *xten=(char *)malloc(header.nsymbt);
  fread(xten,header.nsymbt,1,pFile);


  cout<<"NX :"<<header.nx<<endl;
  cout<<"Ny :"<<header.ny<<endl;
  cout<<"NZ :"<<header.nz<<endl;
  cout<<"xorigin :"<<header.xorigin<<endl;
  cout<<"yorigin :"<<header.yorigin<<endl;
  cout<<"zorigin :"<<header.zorigin<<endl;
  cout<<"mode    :"<<header.mode   <<endl;
  cout<<"mapc    :"<<header.mapc   <<endl;
  cout<<"mapr    :"<<header.mapr   <<endl;
  cout<<"maps    :"<<header.maps   <<endl;
  cout<<"xten    :"<<xten          <<endl;
  int ord = 0;
  char * m2=(char *)&header;
  if (m2[0]==0 && m2[1]==0) ord=1; // int's are typically 32-bit, this one should be nx.  
  cout<<"file endianness = "<<ord<<", where 1 = big-endian, 0 = little-endian"<<endl;
  

  long x = 0x34333231;
  char *y = (char *) &x;
  int mord;
  if(strncmp(y,"1234",4)){
    cout<<"system is Big Endian"<<endl;
    mord = 1; }
  else {
    cout<< "system is little Endian"<<endl;
    mord = 0; }
  if  (mord ^ ord) 
  {
    cout<<"REVERSING BYTE ORDER"<<endl;
    for (int i=0; i<56*4; i+=4) 
    { 
      char w=m2[i]; m2[i]=m2[i+3]; m2[i+3]=w; 
      w=m2[i+1]; m2[i+1]=m2[i+2]; m2[i+2]=w;
    }
  }
  //header = mrcH(m2); 



  cout<<"NX :"<<header.nx<<endl;
  cout<<"Ny :"<<header.ny<<endl;
  cout<<"NZ :"<<header.nz<<endl;
  cout<<"xorigin :"<<header.xorigin<<endl;
  cout<<"yorigin :"<<header.yorigin<<endl;
  cout<<"zorigin :"<<header.zorigin<<endl;
  cout<<"mode    :"<<header.mode   <<endl;
  cout<<"mapc    :"<<header.mapc   <<endl;
  cout<<"mapr    :"<<header.mapr   <<endl;
  cout<<"maps    :"<<header.maps   <<endl;
  cout<<"xten    :"<<xten          <<endl;


  //unsigned char* cdata[1][1][1];//[header.nx][header.ny][header.nz];
  unsigned char* cdata[header.nx][header.ny][header.nz];
  //fseek (pFile , 0 , SEEK_SET);
      //fread (&cdata, 1            , 1 , pFile);
      //cout<<"results of fread: "<<fread (&cdata, 1            ,1  , pFile)<<endl;;
      cout<<"results of fread: "<<fread (&cdata, 1            , header.nx* header.ny* header.nz , pFile)<<endl;;
      cout<< "cdata :"<<cdata<<endl;
      cout<< "*cdata :"<<*cdata<<endl;
  for (int k = 0; k < header.nz; k++ ) for ( int j = 0 ; j  < header.ny; j++) for (int i = 0 ; i < header.nx; i++)  {
      cout<<" i : "<<i<<endl;
      cout<<" j : "<<j<<endl;
      cout<<" k : "<<k<<endl;
      cout<<(cdata[i][j][k])<<endl;
      //fread (&cdata, 1            ,1, pFile);
      //cout<< "cdata :"<<cdata<<endl;


  }

  //result = fread (buffer,1,1024 ,pFile); 
  //was 
  //result = fread (buffer,1,lSize,pFile); 
  //cout << buffer<<std::endl;
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

  /* the whole file is now loaded in the memory buffer. */

  // terminate
  fclose (pFile);
  free (buffer);
  return 0;
}
