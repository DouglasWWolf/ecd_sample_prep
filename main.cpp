//=================================================================================================
// ecd_sample_prep  
//
// Author: Doug Wolf
//
// Command line options:
//
//   -config <filename>    : specifies the name of a configuration file
//
//   -trace <cell_number>  : instead of creating an output file, traces a cell in an existing 
//                           file.
//                        
//=================================================================================================

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <fcntl.h>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <fstream>
#include "config_file.h"

using namespace std;

void     execute(const char** argv);
void     loadFragments();
void     loadDistribution();
uint32_t findLongestSequence();
uint32_t verifyDistributionIsValid();
void     writeOutputFile(uint32_t frameGroupCount);
void     parseCommandLine(const char** argv);
void     trace(uint32_t cellNumber);
void     readConfigurationFile(string filename);

// Define a convenient type to encapsulate a vector of strings
typedef vector<string> strvec_t;

// Contains nucleic acid fragement definitions
map<string, vector<int>> fragment;

// This list defines each fragment distribution in the distribution definitions file
struct distribution_t
{
    int             first, last, step;
    vector<uint8_t> cellValue;
};
vector<distribution_t> distributionList;

//=================================================================================================
// Command line options
//=================================================================================================
struct cmdline_t
{
    bool     trace;
    uint32_t cellNumber;
    string   config;

} cmdLine;
//=================================================================================================


//=================================================================================================
// Variable names in this structure should exactly match the configuration file
//=================================================================================================
struct config_t
{
    uint32_t cells_per_frame;
    uint64_t contig_size;
    uint8_t  diagnostic_constant;
    uint32_t diagnostic_frames;
    uint32_t data_frames;
    uint8_t  quiescent;
    string   fragment_file;
    string   distribution_file;
    string   output_file;

} config;
//=================================================================================================



//=================================================================================================
// main() - Execution starts here.
//=================================================================================================
int main(int argc, const char** argv)
{
    try
    {
        execute(argv);
    }
    catch(const std::exception& e)
    {
        cerr << e.what() << '\n';
    }
}
//=================================================================================================



//=================================================================================================
// throwRuntime() - Throws a runtime exception
//=================================================================================================
static void throwRuntime(const char* fmt, ...)
{
    char buffer[1024];
    va_list ap;
    va_start(ap, fmt);
    vsprintf(buffer, fmt, ap);
    va_end(ap);

    throw runtime_error(buffer);
}
//=================================================================================================

//=================================================================================================
// parseCommandLine() - Parse the command line parameters into the -cmdLine structure
//=================================================================================================
void parseCommandLine(const char** argv)
{
    int i=0;

    // Loop through each command line parameter
    while (argv[++i])
    {
        // Fetch this parameter
        string token = argv[i];

        // Handle the "-trace" command line switch
        if (token == "-trace")
        {
            cmdLine.trace = true;
            if (argv[i+1])
                cmdLine.cellNumber = atoi(argv[++i]);                
            else
                throwRuntime("Missing parameter on -trace");                
            continue;
        }

        // Handle the "-config" command line switch
        if (token == "-config")
        {
            if (argv[i+1])
                cmdLine.config = argv[++i];
            else
                throwRuntime("Missing parameter on -config");
            continue;
        }

        printf("Illegal command line parameter '%s'\n", token.c_str());
        exit(1);
    }
}
//=================================================================================================




//=================================================================================================
// execute() - Top level code for program logic
//=================================================================================================
void execute(const char** argv)
{
    // Ensure that comma-separators get printed for numbers
    setlocale(LC_ALL, "");

    // Parse the command line
    parseCommandLine(argv);

    // Fetch the configuration values from the file and populate the global "config" structure
    readConfigurationFile(cmdLine.config);

    // If we're supposed to trace a single cell, make it so
    if (cmdLine.trace)
    {
        trace(cmdLine.cellNumber);
        exit(0);
    }

    // Load the fragment definitions
    loadFragments();

    // Load the fragment sequence distribution definitions
    loadDistribution();

    // Find out how many frame groups we need to write to the output file
    uint32_t frameGroupCount = verifyDistributionIsValid();

    // Write the output file
    writeOutputFile(frameGroupCount);
}
//=================================================================================================



//=================================================================================================
// getNextCommaSeparatedToken() - Fetches the next comma separated token from a line of text
//
// Passed:  p (by REFERENCE!) = A pointer to the text being parsed
//          token             = A pointer to where the extracted token should be stored
//
// Returns: true if a token was extracted, false if no more tokens available on the line
//
// Note: This routine will accomodate lines containing optional carriage-returns
//       Line is assumed to contain no linefeeds
//=================================================================================================
bool getNextCommaSeparatedToken(const char*& p, char* token)
{
    // Clear the caller's 'token' field in case we can't find a token
    *token = 0;

    // Skip over white-space
    while (*p == 32 || *p == 9) ++p;

    // If we've hit the end of the input line, tell the caller
    if (*p == 0 || *p == 10 || *p == 13) return false;

    // Extract the token into the buffer
    while (!(*p == 32 || *p == 9 || *p == 10 || *p == 13 || *p == 0 || *p == ',')) *token++ = *p++;

    // Nul-terminate the extracted token
    *token = 0;

    // Skip over any trailing whitespace
    while (*p == 32 || *p == 9) ++p;

    // If there is a comma, skip over it
    if (*p == ',') ++p;

    // Tell the caller that they have extracted a token
    return true;
}
//=================================================================================================


//=================================================================================================
// getNextCommaSeparatedInt() - Similar to "getNextCommaSeparatedToken()", but converts the token
//                              to an integer
//=================================================================================================
bool getNextCommaSeparatedInt(const char*& p, int* pValue)
{
    char token[1000];

    // Fetch the next token
    bool status = getNextCommaSeparatedToken(p, token);    

    // Convert that token to an integer
    *pValue = atoi(token);

    // And tell the caller whether or not there was a token available
    return status;
}
//=================================================================================================


//=================================================================================================
// loadFragments() - Load fragment definitions into RAM
//
// On Exit: the global "fragment" object contains fragment definitions
//=================================================================================================
void loadFragments()
{
    char fragmentName[1000], buffer[1000];
    vector<int> v;
    string line;

    // Fetch the filename of the fragment definiton file
    const char* filename = config.fragment_file.c_str();

    // Open the input file
    ifstream file(filename);   

    // If we can't open the input file, complain
    if (!file.is_open()) throwRuntime("%s not found", filename);

    // Loop through each line of the input file
    while (getline(file, line))
    {
        // Get a pointer to the line of text we just read
        const char* p = line.c_str();

        // Skip over whitespace
        while (*p == 32 || *p == 9) ++p;

        // If the line is blank, skip it
        if (*p == 0 || *p == 10 || *p == 13) continue;

        // Any line starting with '#' is a comment
        if (*p == '#') continue;

        // Any line starting with '//' is a comment
        if (p[0] == '/' && p[1] == '/') continue;

        // Clear the fragment value vector
        v.clear();

        // Fetch the fragment name
        getNextCommaSeparatedToken(p, fragmentName);

        // If the fragment name is blank, skip this line
        if (fragmentName[0] == 0) continue;

        // Fetch every integer value after the name
        while (getNextCommaSeparatedToken(p, buffer))
        {
            v.push_back(atoi(buffer));
        }

        // Save this fragment data into our global variable
        fragment[fragmentName] = v;
    }
}
//=================================================================================================




//=================================================================================================
// dumpDistributionList() - Displays the distribution list for debugging purposes
//=================================================================================================
void dumpDistributionList()
{
    for (auto& r : distributionList)
    {
        printf("%i : %i : %i  *** ", r.first, r.last, r.step);
        for (auto i : r.cellValue) printf("%d  ", i);
        printf("\n");

    }
}
//=================================================================================================



//=================================================================================================
// loadDistribution() - Loads the fragment distribution definitions into RAM
//
// On Exit: the global "distribuitionList" object contains distribution definitions
//=================================================================================================
void loadDistribution()
{
    char fragmentName[1000];
    distribution_t distRecord;
    string line;

    // Get a handy reference to the vector of cell values in a distribution record
    auto& drcv = distRecord.cellValue;

    // Fetch the filename of the fragment distribiution definiton file
    const char* filename = config.distribution_file.c_str();

    // Open the input file
    ifstream file(filename);   

    // If we can't open the input file, complain
    if (!file.is_open()) throwRuntime("%s not found", filename);

     // Loop through each line of the input file
    while (getline(file, line))
    {
        // Get a pointer to the line of text we just read
        const char* p = line.c_str();

        // Skip over whitespace
        while (*p == 32 || *p == 9) ++p;

        // If the line is blank, skip it
        if (*p == 0 || *p == 10 || *p == 13) continue;

        // Any line starting with '#' is a comment
        if (*p == '#') continue;

        // Any line starting with '//' is a comment
        if (p[0] == '/' && p[1] == '/') continue;

        // Look for the '$' delimeter that begins a list of fragment IDs
        const char* delimeter = strchr((char*)p, '$');

        // If that delimeter doesn't exist, this isn't a valid distribution definition
        if (delimeter == nullptr) continue;

        // Replace the '$' with a carriage return to fake an "end of line"
        *(char*)delimeter++ = 13;

        // Just in case the user added a comma after the '$', consume it 
        while (*delimeter == 32 || *delimeter == 9) ++delimeter;
        if (*delimeter == ',') ++delimeter;

        // Get the first cell number, the last cell number, and the step-size
        getNextCommaSeparatedInt(p, &distRecord.first);
        getNextCommaSeparatedInt(p, &distRecord.last );
        getNextCommaSeparatedInt(p, &distRecord.step );

        // Ensure that the first cell number in the distribution is valid
        if (distRecord.first < 1 || distRecord.first > config.cells_per_frame)
        {
            throwRuntime("Invalid cell number %i", distRecord.first);
        }

        // If no "last cell" was specified, this distribution is just for the first cell
        if (distRecord.last == 0) distRecord.last = distRecord.first;

        // If no 'step' is specified, we're defining every cell from 'first' to 'last'
        if (distRecord.step == 0) distRecord.step = 1;

        // Clear the vector that will hold fragment data values
        drcv.clear();

        // Point to the comma separated fragement ids that come after the '$' delimeter
        p = delimeter;

        // Loop through every fragment name in the comma separated list...
        while (getNextCommaSeparatedToken(p, fragmentName))
        {
            // If we don't recognize this fragment name, complain
            if (fragment.find(fragmentName) == fragment.end())
            {
                throwRuntime("Undefined fragment name '%s'", fragmentName);
            }

            // Get a reference to the cell values for this fragment
            auto& fragcv = fragment[fragmentName];

            // Append the cell values for this fragment to the distribution record
            drcv.insert(drcv.end(), fragcv.begin(), fragcv.end());
        }

        // And add this distribution record to the distribution list
        distributionList.push_back(distRecord);
    }
}
//=================================================================================================


//=================================================================================================
// findLongestSequence() - Finds and returns the number of frames requires by the longest sequence
//                         in the distributionList
//=================================================================================================
uint32_t findLongestSequence()
{
    uint32_t longestLength = 0;

    // Loop through every record in the distribution list and keep
    // track of the length of the longest sequence of fragments we find    
    for (auto& distRec : distributionList)
    {
        // Keep track of the length of the longest sequence of fragments we find
        if (distRec.cellValue.size() > longestLength) longestLength = distRec.cellValue.size();
    };

    // Hand the caller the length of the longest sequence of fragments
    return longestLength;
}
//=================================================================================================


//=================================================================================================
// verifyDistributionIsValid() - Checks to make sure that number of frame groups implied by the
//                               longest fragement sequence will fit into the contiguous buffer.
//
// Returns: The number of frame groups that will be written to the output file
//=================================================================================================
uint32_t verifyDistributionIsValid()
{
    // What's the maximum number of frames that will fit into the contig buffer?
    uint32_t maxFrames = config.contig_size / config.cells_per_frame;

    // What is the maximum number of frames required by any fragment sequence?
    uint32_t longestSequence = findLongestSequence();

    // A "frame group" is a set of diagnostic frames followed by a set of data frames.
    uint32_t frameGroupLength = config.diagnostic_frames + config.data_frames;

    // How many frames groups are required to express our longest sequence?
    uint32_t frameGroupCount = longestSequence / config.data_frames + 1;

    // How many data frames are in 'frameGroupCount' frameg groups?
    uint32_t totalReqdFrames = frameGroupCount * frameGroupLength;

    // How many bytes will that number of frames occupy in the contiguous buffer?
    uint64_t totalContigReqd = (uint64_t)totalReqdFrames * (uint64_t)config.cells_per_frame;

    // Tell the user basic statistics about this run
    printf("%'16u Frames in the longest fragment sequence\n", longestSequence);
    printf("%'16u Frames in a frame group\n", frameGroupLength);
    printf("%'16u Frame group(s) required\n", frameGroupCount);
    printf("%'16u Frames will fit into the contig buffer\n", maxFrames);
    printf("%'16u Frames required in total\n", totalReqdFrames);
    printf("%'16lu Bytes required in total\n", totalContigReqd);

    // If the longest fragment sequence is too long to fit into the contiguous buffer,
    // complain and drop dead
    if (totalReqdFrames > maxFrames)
    {
        printf("\nThe specified fragment distribution won't fit into the contiguous buffer!\n");
        exit(1);
    }

    // Tell the caller how many frame groups we're going to output
    return frameGroupCount;
}
//=================================================================================================


//=================================================================================================
// buildDataFrame() - Uses the fragment-sequence distribution list to create a data frame
//=================================================================================================
void buildDataFrame(uint8_t* frame, uint32_t frameNumber)
{
    // Every cell in the frame starts out quiescient
    memset(frame, config.quiescent, config.cells_per_frame);

    // Loop through every distribution record in the distribution list
    for (auto& dr : distributionList)
    {
        // If this fragment sequence contains a value for this frame number...
        if (frameNumber < dr.cellValue.size())
        {
            // Populate the appropriate cells with the data value for this frame
            for (uint32_t cellNumber = dr.first-1; cellNumber < dr.last; cellNumber += dr.step)
            {
                frame[cellNumber] = dr.cellValue[frameNumber];
            }
        }
    }
}
//=================================================================================================

//=================================================================================================
// writeOutputFile() - Creates the output file
//=================================================================================================
void writeOutputFile(uint32_t frameGroupCount)
{
    uint32_t i, frameNumber = 0;

    // Fetch the name of the file we're going to create
    const char* filename = config.output_file.c_str();
   
    // Open the file we're going to write, and complain if we can't
    FILE* ofile = fopen(filename, "w");
    if (ofile == nullptr) throwRuntime("Can't create %s", filename);

    // Allocate sufficient RAM to contain an entire data frame
    unique_ptr<uint8_t> framePtr(new uint8_t[config.cells_per_frame]);

    // Get a pointer to the frame data
    uint8_t* frame = framePtr.get();

    // Loop through each frame group
    for (int32_t frameGroup = 0; frameGroup < frameGroupCount; ++frameGroup)
    {
        // Build the diagnostic frame
        memset(frame, config.diagnostic_constant, config.cells_per_frame);

        // Write the correct number of diagnostic frames to the output file
        for (i=0; i<config.diagnostic_frames; ++i)
        {
            fwrite(frame, 1, config.cells_per_frame, ofile);
        }

        // Write the correct number of data frames to the output file
        for (i=0; i<config.data_frames; ++i)
        {
            buildDataFrame(frame, frameNumber++);
            fwrite(frame, 1, config.cells_per_frame, ofile);
        }
    }

    // We're done with the output file
    fclose(ofile);
}
//=================================================================================================


//=================================================================================================
// trace() - Displays the value of a single cell for every frame in the output file
//=================================================================================================
void trace(uint32_t cellNumber)
{
    bool first = true;

    // Fetch the name of the file we're going to open
    const char* filename = config.output_file.c_str();

    // Open the file we're going to read, and complain if we can't
    FILE* ifile = fopen(filename, "r");
    if (ifile == nullptr) throwRuntime("Can't create %s", filename);

    // Allocate sufficient RAM to contain an entire data frame
    unique_ptr<uint8_t> framePtr(new uint8_t[config.cells_per_frame]);

    // Get a pointer to the frame data
    uint8_t* frame = framePtr.get();

    // Loop through each frame of the file...
    while (fread(frame, 1, config.cells_per_frame, ifile) == config.cells_per_frame)
    {
        // If this isn't the first value we've output, print a comma separator
        if (!first) printf(", ");
        first = false;
        
        // And display the value of the cell number that was specified on the command line
        printf("%d", frame[cellNumber]);
    }
    
    // Terminate the line of text in the output
    printf("\n");
}
//=================================================================================================


//=================================================================================================
// readConfigurationFile() - Reads in the configuration file and populates the global "config"
//                           structure.
//=================================================================================================
void readConfigurationFile(string filename)
{
    CConfigFile cf;

    // Declare a default filename
    const char* cfilename = "ecd_sample_prep.conf";

    // If the filename passed by the caller isn't blank, that's our filename
    if (!filename.empty()) cfilename = filename.c_str();

    // Read and parse the configuration file and complain if we can't
    if (!cf.read(cfilename, false)) throwRuntime("Can't read %s", cfilename);

    // Fetch each configuration
    cf.get("cells_per_frame",     &config.cells_per_frame    );
    cf.get("contig_size",         &config.contig_size        );
    cf.get("diagnostic_frames",   &config.diagnostic_frames  );
    cf.get("data_frames",         &config.data_frames        );
    cf.get("diagnostic_constant", &config.diagnostic_constant);
    cf.get("quiescent",           &config.quiescent          );
    cf.get("fragment_file",       &config.fragment_file      );
    cf.get("distribution_file",   &config.distribution_file  );
    cf.get("output_file",         &config.output_file        );
}
//=================================================================================================
