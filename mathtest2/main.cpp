//
//  main.cpp
//  mathtests
//
//  Created by Harvard Virgil Humphrey on 2019/7/22.
//  Copyright © 2019 Harvard Virgil Humphrey. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include "arrayfunctions.h"
#include "matrixfunctions.h"
#include "mathfunctions.h"
#include "rosless.hpp"
//#include <eigen3/Eigen/Dense>

int main() {
    double vec[12] = {-0.05,0.05,0.4,0.5,-0.6,0.9,0.7,-0.75,0.1,-0.4,-1.2,-0.4};
    double ang[6] = {2.8928,2.8853,-0.6270,-0.8390,-2.3967,0.5347};
    
    double** state;
    zerosdouble(state,6,2);
    double* angs = new double[6];
    int counter = 0;
    for(int i=0;i<6;i++)    {
        angs[i] = ang[i];
        for(int j=0;j<2;j++)    {
            state[i][j] = vec[counter];
            counter++;
        }
    }
    
    barrier test(6);
    test.populatestructs(state,angs,6);
    test.update_sparse();
    //printdoublemat(test.Aprox,6,6);
    test.calculate_u();
    test.unicycle_dynamics();
    //printdoublemat(test.u_data,6,2);
    printdoublearray(test.input_data.v,6);
    printdoublearray(test.input_data.w,6);


    return 0;
    
}
