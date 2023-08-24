#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <unordered_set>
#include <map>
using namespace std;

// compile with g++ -std=c++11
//
// ./norm_bed reads.bed file.bed colNum [ -len ] > output.bed
//
//

// Function to normalize read counts by region length and library size 
int normBed( std::string readsName, std::string inBedName, int colNum, bool lengthNorm ) {
    std::ifstream readsFile( readsName );
    std::ifstream inBedFile( inBedName );

    // Count lines in readsFile
    std::string readsLine;
    double readCount = 0;

    while ( std::getline(readsFile, readsLine) ) {
        readCount += 1;
    }

    // Parse inBedFile
    std::string bedLine;
    std::string delim = "\t";
    size_t      pos   = 0;

    while ( std::getline(inBedFile, bedLine) ) {
        std::vector<std::string> cols;

        // Create vector containing all columns for the line
        while ( (pos = bedLine.find(delim)) != std::string::npos ) {
            std::string col = bedLine.substr( 0, pos );
            cols.push_back( col );
            bedLine.erase( 0, pos + delim.length() );
        }
                  
        // Add last column to the vector
        cols.push_back( bedLine );

        // Assign variables for each column
        std::string chrom = cols[0];
        int         start = std::atoi( cols[1].c_str() );
        int         end   = std::atoi( cols[2].c_str() );
        double      count = std::atof( cols[ colNum - 1 ].c_str() );

        // Normalize by interval length
        if ( lengthNorm == true ) {
            count = count / (end - start);            
        }

        // Normalize by library size
        count = count / (readCount / 1000000);

        // Add normalized counts back to column vector
        cols[ colNum - 1 ] = std::to_string(count);

        // Output columns
        for ( std::vector<std::string>::iterator it = cols.begin(); it != cols.end() - 1; ++it) {
            std::cout << (*it) << "\t";
        }

        int lastCol = cols.size() - 1;
        std::cout << cols[ lastCol ] << endl;
    }
}



// main function
int main( int argc, char *argv[] ) {
    std::string readsName = argv[1];
    std::string inBedName = argv[2];
    std::string inColNum = argv[3];

    int colNum = std::atoi( inColNum.c_str() );

    bool lengthNorm = false;

    if ( argc == 5 ) {
        std::string lengthArg = argv[4];

        if ( lengthArg == "-len" ) {
            lengthNorm = true;
        }
    }

    normBed( readsName, inBedName, colNum, lengthNorm );

    return 0;
}



