//
//  utility.h
//  Master Thesis
//
//  Created by Magnus Mencke on 30/03/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef utility_hpp
#define utility_hpp

#include <fstream>
#include <iostream>
#include <ql/math/matrix.hpp>

#include <ql/models/model.hpp>

using namespace QuantLib;

inline void writeCsv(const Matrix data, std::string file_name) {
    std::ofstream myfile;
    myfile.open (file_name);
    
    Size n =data.rows();
    Size m = data.columns();
    
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            myfile << data[i][j];
            if(j+1<m) {
                myfile << ", ";
            }
        }
        myfile << "\n";
    }
    
    myfile.close();
}


inline Matrix readCsv(std::string file_name) {
    std::ifstream myfile;
    myfile.open(file_name);
    
    std::string myline,myentry,temp;
    
    std::vector<std::string> row;
    
    std::vector<std::vector<std::string>> matrix;
    
    if(myfile.good()) {
        // read an entire row and
        // store it in a string variable 'line'
        while(getline(myfile,myline)) {
            
            // used for breaking words
            std::stringstream ss(myline);
            
            row.clear();
            
            while (getline(ss, myentry, ',')) {
                if(myentry=="\r") {
                    //don't add element if Carriage return
                    //(might need to make this more robust if I run into issues again)
                } else {
                    row.push_back(myentry);
                }
            }
            matrix.push_back(row);
        }
    }
    
    
    size_t n=matrix.size();
    size_t m=matrix[0].size();
    Matrix data(n,m);
    
    for(int i=0; i<n;i++) {
        for(int j=0;j<m;j++) {
            data[i][j]=std::stold(matrix[i][j]);//convert to number
        }
    }
    return(data);
}


inline void newSection(std::string section_title) {
    
    std::string dash(section_title.length(), '-');
    
    std::cout << std::endl << "|-" << dash << "-|" << std::endl;
    std::cout << "| " << section_title << " |" << std::endl;
    std::cout << "|-" << dash << "-|" << std::endl << std::endl;
}


inline Rate discountRate(AffineModel &the_model, Time start, Time end) {
    Real result;
    
    
    if(end-start>0) {
        result = -std::log(the_model.discount(end)/the_model.discount(start))/(end-start);
    } else {
        result = -std::log(the_model.discount(end))/(end);
    }
    
    return result;
}


inline void printMatrix(const Matrix mat) {
    for(int i=0;i<mat.rows();i++) {
        for(int j=0;j<mat.columns();j++) {
            std::cout << mat[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


#endif /* utility_hpp */
