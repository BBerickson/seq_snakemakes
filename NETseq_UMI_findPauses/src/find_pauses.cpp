#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <unordered_set>
#include <thread>
#include <chrono>
#include <sys/stat.h>
using namespace std;

// compile with g++ -std=c++11 find_pauses.cpp -o find_pauses



// Helper to create vector containing all columns for the line
std::vector<std::string> createDelimVec(std::string lineIn, std::string delim) {
    std::vector<std::string> cols;

    size_t pos = 0;

    while ((pos = lineIn.find(delim)) != std::string::npos) {
        std::string col = lineIn.substr(0, pos);
        cols.push_back(col);
        lineIn.erase(0, pos + delim.length());         
    }

    cols.push_back(lineIn);  // add last column to vector

    return cols;
}





// Split bedGraph into single basepair intervals 
int splitBedGraph(std::string bedGraphName, std::string bedGraphSplit) {

    std::ifstream bedGraph(bedGraphName);
    ofstream outputFile(bedGraphSplit.c_str());
    
    std::string line;
 
    while (std::getline(bedGraph, line)) {
        std::string lineCol = line;

        std::vector<std::string> cols = createDelimVec(lineCol, "\t");

        // Assign variables for each column
        std::string chrom = cols[0];
        int         start = std::atoi(cols[1].c_str());
        int         end   = std::atoi(cols[2].c_str());
        std::string count = cols[3];

        int length = end - start;

        // Divide bedGraph into single bp intervals 
        for (int i = 0; i < length; i++) {
            int startNew = start + i;
            int endNew   = startNew + 1;
        
            std::string startOut = std::to_string(startNew);
            std::string endOut   = std::to_string(endNew);

            std::string newLine = 
                chrom    + "\t" + 
                startOut + "\t" + 
                endOut   + "\t" + 
                count    + "\n";
            
            outputFile << newLine;
        }
    }

    bedGraph.close();

    return 0;
}





// Generate windows around each bedGraph interval
// windows that fall outside intervals in geneRegion are adjusted so they
// are entirely within a region and of length winSize
// windows are only created for intervals with signal > countLim
// Output must be sorted based on window name
int runBedtools(std::string bedGraphSplit, std::string geneRegion, std::string chromSizes, 
                int winSize, std::string countLim, std::string threads, std::string bedGraphWins) {

    // bedGraphSplit should follow standard bed format and have the following
    // columns:
    // 1. chrom
    // 2. start
    // 3. end
    // 4. count

    // geneRegion should follow standard bed format and have the following
    // columns:
    // 1. chrom
    // 2. start
    // 3. end
    // 4. name
    // 5. score
    // 6. strand

    // chromSizes should be tab delimited and have the following columns:
    // 1. chrom
    // 2. chrom length

    int leftLen  = winSize / 2;
    int rightLen = leftLen - 1;

    std::string winStr   = std::to_string(winSize);
    std::string leftStr  = std::to_string(leftLen);
    std::string rightStr = std::to_string(rightLen);

    std::string bashCommands =
        // get bedGraph signal that intersects with geneRegion
        // this produces a bed file with the following columns:
        // 1.  gene chrom
        // 2.  gene start
        // 3.  gene end
        // 4.  gene name
        // 5.  gene score
        // 6.  gene strand
        // 7.  bg chrom
        // 8.  bg start
        // 9.  bg end
        // 10. bg count
        // 11. overlap length
        "cat " + geneRegion +
        "    | sort -k1,1 -k2,2n " +
        "    | bedtools intersect -sorted -wo -a - -b " + bedGraphSplit +

        // filter intervals and format columns
        // filter for candidate pause sites by filtering for intervals where
        // signal > countLim
        // create window name by combining the candidate pause site coordinates,
        // gene coordinates, and gene name
        // this produces a bed file with the following columns:
        // 1. pause chrom
        // 2. pause start
        // 3. pause end
        // 4. window name
        // 5. pause score
        // 6. pause strand
        // 7. pause start
        // 8. pause end
        // 9. pause count
        "    | awk -v OFS=\"\\t\" '($10 >= " + countLim + ") { " +
        "        chr   = $1; " +
        "        nm    = $4; " +
        "        score = $5; " +
        "        strnd = $6; " +
        "        strt  = $8; " +
        "        end   = $9; " +
        "        count = $10; " +
        "        nm    = chr\":\"strt\"-\"end\"*\"nm; " +

        "        print chr, strt, end, nm, score, strnd, strt, end, count " +
        "    }' " +

        // create windows by expanding candidate pause coordinates based on
        // winSize
        "    | bedtools slop -i - -g " + chromSizes + " -l " + leftStr + " -r " + rightStr +
        "    | sort -S1G --parallel=" + threads + " -k1,1 -k2,2n -k3,3n -u " +

        // get bedGraph signal that intersects with each window
        // the output must be sorted based on window name (column 4)
        // this produces a bed file with the following columns:
        // 1.  win chrom
        // 2.  win start
        // 3.  win end
        // 4.  win name
        // 5.  win score
        // 6.  win strand
        // 7.  pause start
        // 8.  pause end
        // 9.  pause count
        // 10. bg chrom
        // 11. bg start
        // 12. bg end
        // 13. bg count
        // 14. overlap length
        "    | bedtools intersect -sorted -a - -b " + bedGraphSplit + " -wo " +
        "    | sort -S1G --parallel=" + threads + " -k4,4 -k1,1 -k2,2n -k3,3n " +
        "    > " + bedGraphWins;

    system(bashCommands.c_str());

    return 0;
}





// Check pause coordinates and pause count to catch malformed input data
int checkCountMatch(int pauseStart, int bedStart, int pauseCount, int count) {
    if (pauseStart == bedStart && pauseCount != count) {
        std::cerr << "ERROR: malformed input data, bedGraph count "
                  << "does not match pauseCount when coordinates "
                  << "match"
                  << endl;
        exit(1);
    }

    return 0;
}

int checkCountLim(int pauseCount, int countLim) {
    if (pauseCount < countLim) {
        std::cerr << "ERROR: malformed input data, pauseCount is "
                  << "less than countLim, only candidate pauses "
                  << "with counts >= countLim should be included "
                  << "in input data"
                  << endl;
        exit(1);
    }

    return 0;
}

int pauseErrorMsg(std::string str1, std::string str2) {
    std::cerr << "ERROR: malformed input data " + str1
              << " does not equal " + str2
              << endl;
    exit(1);
}

int checkPauseCoords(int oldPauseCount, int pauseCount, int oldPauseStart,
                     int pauseStart, int oldPauseEnd, int pauseEnd) {

    if (oldPauseCount != pauseCount) {
        pauseErrorMsg("oldPauseCount", "pauseCount");
    }

    if (oldPauseStart != pauseStart) {
        pauseErrorMsg("oldPauseStart", "pauseStart");
    }

    if (oldPauseEnd != pauseEnd) {
        pauseErrorMsg("oldPauseEnd", "pauseEnd");
    }

    return 0;
}





// Calculate signal cutoff given the window size and a vector of counts
double calcWinLimit(int winSize, std::vector<double> counts, int stdevLim) {

    // Number of non-zero values
    int countLen = counts.size();

    // Sum counts
    double countTot = 0;

    for (int i = 0; i < countLen; i++) {
        countTot += counts[i];
    }
    
    // Calculate number of zeros to add to counts
    double zeroLen = winSize - countLen;

    // Add zeros for bedGraph intervals that lacked signal
    for (int i = 0; i < zeroLen; i++) {
        counts.push_back(0);
    }

    // Calculate the mean and standard deviation for each window
    double meanCount = countTot / winSize;
    double sumVar = 0;

    for (int i = 0; i < counts.size(); i++) {
        double var = counts[i] - meanCount; 
        sumVar += pow(var, 2);
    }
    
    double meanVar = sumVar / (winSize - 1);
    double stdev   = sqrt(meanVar);
    double limit   = meanCount + (stdev * stdevLim);

    return limit;
}





// Print pause coordinates and metrics
std::unordered_set<std::string> writePause(std::string chrom, int start, int end, std::string name,
                                           std::string strand, int round, int pauseCount, int totCount,
                                           std::unordered_set<std::string> pauseSet, double limit) {
    // Pause stats
    std::string roundStr      = std::to_string(round);
    std::string pauseCountStr = std::to_string(pauseCount);
    std::string totCountStr   = std::to_string(totCount);
    
    std::string pauseStats = roundStr + ":" + pauseCountStr + "," + totCountStr;

    // Write pause coords to stdout
    std::cout << chrom      << "\t"
              << start      << "\t" 
              << end        << "\t"
              << name       << "\t" 
              << pauseStats << "\t"
              << strand     << endl;

    // Add pause coords to pauseSet
    std::string startStr = std::to_string(start);
    std::string endStr   = std::to_string(end);

    std::string pauseCoords = chrom + ":" + startStr + "-" + endStr + strand;

    pauseSet.insert(pauseCoords);

    return pauseSet;
}





// Identify pause sites
int findPauses(std::string bedGraphWins, int countLim, int winLim, int stdevLim) {

    // To identify pauses:
    // 1. The pause counts must be > countLim
    // 2. The pause counts must be stdevLim standard deviations greater than
    //    the mean window counts
    // 3. The window counts (minus pause counts) must be > winLim
 
    std::unordered_set<std::string> pauseSet;  // set to store all pause sites
    int pauseNum = -1;
    int round    = 0;

    while (pauseNum != 0) {
        round   += 1;
        pauseNum = 0;

        std::unordered_set<std::string> roundPauseSet;  // set to store pauses sites from the round

        std::ifstream       Wins(bedGraphWins);
        std::string         oldName    = "";
        double              countTot   = 0;
        int                 winSize    = 0;
        std::string         oldWinChrom;
        int                 oldPauseStart;
        int                 oldPauseEnd;
        int                 oldPauseCount;
        std::string         oldStrand;
        std::vector<double> counts;
        std::string         line;

        // Read windows line by line
        while (std::getline(Wins, line)) {

            std::vector<std::string> cols = createDelimVec(line, "\t");
        
            // Asign variables for each column
            std::string winChrom    = cols[0];
            std::string winStartStr = cols[1];
            std::string winEndStr   = cols[2];
            int         winStart    = std::stoi(winStartStr);
            int         winEnd      = std::stoi(winEndStr);
            std::string name        = cols[3];
            std::string strand      = cols[5];

            std::string pauseStartStr = cols[6];
            std::string pauseEndStr   = cols[7];
            int         pauseStart    = std::stoi(pauseStartStr);
            int         pauseEnd      = std::stoi(pauseEndStr);
            int         pauseCount    = std::stoi(cols[8]); 

            std::string bedStartStr = cols[10];
            std::string bedEndStr   = cols[11];
            int         bedStart    = std::stoi(bedStartStr);
            int         bedEnd      = std::stoi(bedEndStr);

            double      count       = std::stod(cols[12]); 
            double      overlap     = std::stod(cols[13]); 

            winSize = winEnd - winStart;
        
            // CHECK INPUT VALUES
            // Input windows are filtered earlier based on countLim these are
            // redundant checks to catch malformed input data
            checkCountLim(pauseCount, countLim);
            checkCountMatch(pauseStart, bedStart, pauseCount, count);
            
            if (overlap != 1) {
                std::cerr << "ERROR: malformed input data, the overlap between "
                          << "the pause window and bedGraph coordinate "
                          << "(overlap) should always be 1"
                          << endl;
                exit(1);
            }

            // CHECK FOR REDUNDANT PAUSE COORDS
            // Check for pause coords in pauseSet
            // if the pause coordinates have already been identified as a
            // pause site, skip entire window
            std::string pauseCheck = winChrom + ":" + pauseStartStr + "-" + pauseEndStr + strand;

            if (pauseSet.find(pauseCheck) != pauseSet.end()) {
                continue;
            } 

            // CHECK IF COORDS ARE A PAUSE
            // Check for bed coords in pauseSet
            // if the current position has been previously identified as a
            // pause site, exclude from calculations for the window
            std::string bedCheck = winChrom + ":" + bedStartStr + "-" + bedEndStr + strand;

            if (pauseSet.find(bedCheck) != pauseSet.end()) {
                continue;
            }

            // Add up counts that are present in the same window
            if (oldName == "" || name == oldName) {
                counts.push_back(count);
                countTot += count;

                // When the window name is the same the candidate pause
                // coordinates and count should also be the same
                if (name == oldName) {
                    checkPauseCoords(
                        oldPauseCount, pauseCount,
                        oldPauseStart, pauseStart,
                        oldPauseEnd,   pauseEnd
                    );
                }

            // When the window changes calculate stats for the previous window
            } else {
                
                // Check that pauseCount and countTot (minus pauseCount) meet
                // minimum cutoffs
                if (countTot - oldPauseCount >= winLim) {

                    // Idenitify bedGraph intervals where the number of counts is
                    // stdev * stdevLim greater than the mean window counts
                    double limit = calcWinLimit(winSize, counts, stdevLim);
    
                    if (oldPauseCount > limit) {
    
                        pauseNum += 1;

                        roundPauseSet = writePause(
                            oldWinChrom, oldPauseStart, oldPauseEnd, oldName,
                            oldStrand, round, oldPauseCount, int(countTot),
                            roundPauseSet, limit
                        );
                    }
                }

                // Reset everything for the next window
                counts.clear();
                countTot = 0;

                // Begin operating on the new window
                counts.push_back(count); 

                countTot += count;
            }

            oldWinChrom   = winChrom;
            oldPauseStart = pauseStart;
            oldPauseEnd   = pauseEnd;
            oldPauseCount = pauseCount;
            oldName       = name;
            oldStrand     = strand;
        }

        // Check last window
        if (countTot - oldPauseCount >= winLim) {

            double limit = calcWinLimit(winSize, counts, stdevLim);
    
            if (oldPauseCount > limit) {
                pauseNum += 1;

                roundPauseSet = writePause(
                    oldWinChrom, oldPauseStart, oldPauseEnd, oldName,
                    oldStrand, round, oldPauseCount, int(countTot),
                    roundPauseSet, limit
                );
            }
        }

        pauseSet.insert(roundPauseSet.begin(), roundPauseSet.end());

        Wins.close();
    }

    return 0;
}





// Generate string for naming output files 
std::string getName(std::string path) {

    // Create vector of elements in file path
    std::vector<std::string> pathVec = createDelimVec(path, "/");

    std::string lastName = pathVec[pathVec.size() - 1];

    return lastName;
}





// Main
int main(int argc, char *argv[]) {

    // Input arguments
    std::string bedGraph   = argv[1];  // Input bedGraph with raw counts
    std::string geneRegion = argv[2];  // Bed file with region to search
    std::string chromSizes = argv[3];  // File with chromosome sizes
    std::string winSizeIn  = argv[4];  // Window size to use for searching
    std::string countLimIn = argv[5];  // Min counts for a pause site
    std::string winLimIn   = argv[6];  // Min counts for window (not including pause)
    std::string stdevLimIn = argv[7];  // Number of standard deviations above window mean
    std::string threads    = argv[8];  // Number of threads for sorting
    std::string tmpDir     = argv[9];  // Directory to write temporary files

    int countLim = std::atoi(countLimIn.c_str());
    int winLim   = std::atoi(winLimIn.c_str());
    int stdevLim = std::atoi(stdevLimIn.c_str());
    int winSize  = std::atoi(winSizeIn.c_str());

    // Extract bedGraph name  
    std::string bedGraphName = getName(bedGraph);
    
    // Temporary files  
    std::string bedGraphSplit = tmpDir + "/" + bedGraphName + ".split";    
    std::string bedGraphWins  = tmpDir + "/" + bedGraphName + ".wins";

    // Split bedGraph into single base pair intervals 
    splitBedGraph(bedGraph, bedGraphSplit);

    // Generate windows around each bedGraph interval
    if (winSize % 2 != 0) {
        exit(1);
    }

    runBedtools(
        bedGraphSplit,  // BedGraph divided into single bp coordinates
        geneRegion,     // Bed file for gene regions to search for pauses
        chromSizes,     // File containing chromosome sizes
        winSize,        // bp size of windows
        countLimIn,     // Min reads at site to create window
        threads,        // Number of threads to use for sorting
        bedGraphWins    // Output file
    );

    // Identify pause sites 
    findPauses(bedGraphWins, countLim, winLim, stdevLim);

    // Remove temporary files
    std::string rmCommands = "rm " + bedGraphSplit + " " + bedGraphWins; 
    system(rmCommands.c_str()); 

    return 0;
}


