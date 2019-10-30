#include "superimposer.h"

//------------------Superimposer------------------//
//std::pair< std::vector<std::vector<float> > , std::vector<float> >
  // // SCF is only interested in the RMSD at the moment. In the future we can see about returning the full translation-rotation matrix.
float superimposer(std::vector< std::vector<float> > coord0, std::vector< std::vector<float> > coord1, unsigned int natm){

  typedef std::vector< std::vector<float> > matrix;
  typedef std::pair< std::vector< std::vector<float> > , std::vector<float> > return_val;
  //coord0 is the target (probably the input glycine)
  //coord1 is the prediction (alanine backbone minus CB, superpose onto glycine then apply transformation to CB)
  //natm is the number of atoms

  //Given code takes a matrix with the opposite dimensions used, so tranpose 90deg clockwise if necessary
  if (coord0[0].size() == 3){
    coord0 = transpose(coord0);
    coord1 = transpose(coord1);
  }


  matrix mtx;
  mtx.resize(3, std::vector<float>(3, 0.0) );
  std::vector<float> t(3, 0.0);
  std::vector<float> vec(3, 0.0);
  float tolerance = 0.0001;
  float err = 0.0;

  //Error checking
  unsigned int MAXATOM = 2000000; // SCF not sure I really need to worry about an upper limit. Was 200, I increased rather dramatically.
  if (natm > MAXATOM){
    std::cerr << "Error! Increase maxatom" << std::endl;
    throw 0;
  }
  if (natm < 4){
    std::cerr << "Error! Too few atoms (must be 4 or greater)" << std::endl;
    throw 0;
  }
  for (unsigned int i = 0; i < natm; ++i){
    if (std::max(coord0[0][i], coord1[0][i]) > 998){
      std::cerr << "Coordinate error!" << coord0[0][i] << " " << coord1[0][i] << std::endl;
      throw 0;
    }
  }

  //---Center on origin---//
  std::vector<float> x0_center(3, 0.0);
  std::vector<float> x1_center(3, 0.0);
  //For each axis, sum all atom coordinates, divide by natm to get avg d to origin
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      x0_center[j] = x0_center[j] + coord0[j][i];
      x1_center[j] = x1_center[j] + coord1[j][i];
    }
  }
  for (int i = 0; i < 3; ++i){
    x0_center[i] = x0_center[i] / float(natm);
    x1_center[i] = x1_center[i] / float(natm);
  }
  std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" x0_center = "<<x0_center[0]<<", "<<x0_center[1]<<", "<<x0_center[2] <<std::endl;
  std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" x1_center = "<<x1_center[0]<<", "<<x1_center[1]<<", "<<x1_center[2] <<std::endl;
  //Initialize matrices to hold new centered coordinates
  matrix x0;
  x0.resize(3, std::vector<float>(natm, 0.0));
  matrix x1;
  x1.resize(3, std::vector<float>(natm, 0.0));


  //Subtract avg distance from origin from actual coordinates, fill new matrices
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      x0[j][i] = coord0[j][i] - x0_center[j];
      x1[j][i] = coord1[j][i] - x1_center[j];
    }
  }


  //---End Center on origin---//


  matrix aa;
  aa.resize(3, std::vector<float>(3, 0.0) );
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      for (int k = 0; k < 3; ++k){
        aa[j][k] = aa[j][k] + (x1[j][i] * x0[k][i]);
      }
    }
  }


  //Initialize rotation matrix
  matrix rotation;
  rotation.resize(3, std::vector<float>(3, 0));
  rotation[0][0] = 1.0; rotation[1][1] = 1.0; rotation[2][2] = 1.0;


  //---Iterative rotation scheme---//
  int iteration_count = 0;
  bool do51 = true; //"This is a way to deal with those nasty gotos in the FORTRAN code"
  int iflag;
  int ix;
  int iy;
  int iz;
  double sigma;
  double gamma;
  double sg;
  double bb;
  double cc;

  while (true){
    if (do51){
      iflag = 0;
      ix = 0;
    }

    //If the number of iterations exceeds 1000, give up
    ++iteration_count;
    if (iteration_count > 500){ break; }

    iy = ix + 1;
    if (iy == 3){
      iy = 0;
    }
    iz = 3 - ix - iy;

    //What are sigma and gamma?
    sigma = aa[iz][iy] - aa[iy][iz];
    gamma = aa[iy][iy] + aa[iz][iz];
    sg = sqrt( sigma * sigma + gamma * gamma );


    //If root of (sigma^2 + gamma^2) = 0. this is tripped
    if (sg == 0){
      //"Goto 50 in FORTRAN code"
      //If do51 was true at the beginning of this iteration, this will break loop
      ++ix;
      if (iflag == 0){ break; }
      if (ix < 3){ do51 = false; }
      else{ do51 = true; }
      continue;
    }

    sg = 1.0 / sg;
    if (fabs(sigma) < ( tolerance * fabs(gamma)) ){
      //"Goto 50 in FORTRAN code"
      ++ix;
      if (iflag == 0){ break; }
      if (ix < 3){ do51 = false; }
      else{ do51 = true; }

      continue;
    }


    for (int i = 0; i < 3; ++i){
      bb = (gamma * aa[iy][i]) + (sigma * aa[iz][i]);
      cc = (gamma * aa[iz][i]) - (sigma * aa[iy][i]);
      aa[iy][i] = bb * sg;
      aa[iz][i] = cc * sg;
      bb = gamma * rotation[iy][i] + sigma * rotation[iz][i];
      cc = gamma * rotation[iz][i] - sigma * rotation[iy][i];
      rotation[iy][i] = bb * sg;
      rotation[iz][i] = cc * sg;
    }
    iflag = 1;

    //Goto 50 in FORTRAN
    ++ix;
    if (iflag == 0){ break; }
    if (ix < 3){ do51 = false; }
    else{ do51 = true; }

    continue;


  } //End while loop


  //Change x1 coordinates

  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      t[j] = 0.0;
      for (int k = 0; k < 3; ++k){
        //Seemingly a rotation
        t[j] = t[j] + (rotation[j][k] * x1[k][i]);
      }
    }
    for (int j = 0; j < 3; ++j){
      x1[j][i] = t[j];
    }
  }


  //Calculate RMSD
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; i < 3; ++i){
      err = err + pow( (x0[j][i] - x1[j][i]) , 2 );
    }
  }
  err = sqrt(err / natm);

  //"un-center" points
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      coord1[j][i] = x1[j][i] + x0_center[j];
    }
  }

  //Build output matrices
  for (int i = 0; i < 3; ++i){
    t[i] = x0_center[i];
    for (int j = 0; j < 3; ++j){
      mtx[i][j] = rotation[i][j];
      t[i] = t[i] - rotation[i][j] * x1_center[j];
    }
    vec[i] = t[i];
  }


  return_val ret = std::make_pair(mtx, vec);
  //return ret; // SCF is only interested in the RMSD at the moment. In the future we can see about returning the full translation-rotation matrix.
  return err;
}

void superimposer_move(std::vector<float>& x, std::vector< std::vector<float> >& mtx, std::vector<float>& vec){
  std::vector<float> y(3, 0.0);
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      y[i] = y[i] + (x[j] * mtx[i][j]);
    }
  }
  for (int i = 0; i < 3; ++i){
    y[i] = y[i] + vec[i];
  }

  x = y;
  return;
}


std::vector< std::vector<float> > transpose(std::vector< std::vector<float> >& coord0){
  std::vector< std::vector<float> > tmp;
  for (unsigned int i = 0; i < coord0[0].size(); ++i){
    std::vector<float> xyz;
    for (unsigned int j = 0; j < coord0.size(); ++j){
      xyz.push_back(coord0[j][i]);
    }
    tmp.push_back(xyz);
    xyz.clear();
  }
  return tmp;
}


//Matrix output for debugging
void printVec(std::vector< std::vector<float> >& vec){
  for (unsigned int i = 0; i < vec.size(); ++i){
    for (unsigned int j = 0; j < (vec[i]).size(); ++j){
      std::cout << vec[i][j] << " ";
    }
    std::cout << std::endl;
  }
}
