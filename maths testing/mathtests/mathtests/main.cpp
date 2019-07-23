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
#include "rosless.hpp"

int main(int argc, const char * argv[]) {
    double vec[12] = {-0.05,0.05,0.4,0.5,-0.6,0.9,0.7,-0.75,0.1,-0.4,-1.2,-0.4};
    
    double** state;
    zerosdouble(state,6,2);
    int counter = 0;
    for(int i=0;i<6;i++)    {
        for(int j=0;j<2;j++)    {
            state[i][j] = vec[counter];
            counter++;
        }
    }
    
    barrier test(6);
    printdoublemat(state,6,2);
    return 0;
    
}
