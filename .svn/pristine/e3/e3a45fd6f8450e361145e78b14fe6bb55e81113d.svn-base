
Mat33 ConvertAngleAxisToMat33 (Vec4 angleAxis) {
    float x = angleAxis[1];
    float y = angleAxis[2];
    float z = angleAxis[3];
    float s = sin(angleAxis[0]);
    float c = cos(angleAxis[0]);
    float t = 1-cos(angleAxis[0]);

    Mat33 outMat33 ;
    outMat33[0][0] = t*x*x+c;
    outMat33[0][1] = t*x*y-s*z;
    outMat33[0][2] = t*x*z+s*y;
    //outMat33[0][0] = [t*x*x+c, t*x*y-s*z, t*x*z+s*y];
    outMat33[1][0] = t*x*y+s*z;
    outMat33[1][1] = t*y*y+c;
    outMat33[1][2] = t*y*z-s*x;
    outMat33[2][0] = t*x*z-s*y;
    outMat33[2][1] = t*y*z+s*x;
    outMat33[2][2] = t*z*z+c ;
    //outMat33[0] = [[t*x*x+c, t*x*y-s*z, t*x*z+s*y]
    //	[t*x*y+s*z,t*y*y+c,t*y*z-s*x]
    //[t*x*z-s*y,t*y*z+s*x,t*z*z+c]] ;

    return outMat33;	

};
