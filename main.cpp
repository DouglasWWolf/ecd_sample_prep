//=================================================================================================
// ecd_sample_prep  
//
// Author: Doug Wolf
//=================================================================================================

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <fcntl.h>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>

using namespace std;

void     execute();
void     loadFragments();
void     loadDistribution();
uint32_t findLongestSequence();
void     verifyDistributionIsValid();

// Define a convenient type to encapsulate a vector of strings
typedef vector<string> strvec_t;

// Contains nucleic acid fragement definitions
map<string, vector<int>> fragment;

// This list defines each fragment distribution in the distribution definitions file
struct distribution_t {int first, last, step; strvec_t fragments;};
vector<distribution_t> distributionList;

// Defines the sequence of fragments in each cell
strvec_t* distribution;

// Variable names in this structure should exactly match the configuration file
struct config_t
{
    int      cells_per_frame;
    uint64_t contig_size;
    uint32_t diagnostic_frames;
    uint32_t data_frames;
    string   fragment_file;
    string   distribution_file;

} config;

int main()
{

    config.cells_per_frame   = 1024;
    config.contig_size       = 1024 * 5;
    config.diagnostic_frames = 12;
    config.data_frames       = 4525;
    config.fragment_file     = "fragments.csv";
    config.distribution_file = "distribution.csv";

    try
    {
        execute();
    }
    catch(const std::exception& e)
    {
        cerr << e.what() << '\n';
    }

}



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
// execute() - Top level code for program logic
//=================================================================================================
void execute()
{
    // Ensure that comma-separators get printed for numbers
    setlocale(LC_ALL, "");

    loadFragments();

    loadDistribution();

    verifyDistributionIsValid();
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
    if (*p == 0 || *p == 13) return false;

    // Extract the token into the buffer
    while (!(*p == 32 || *p == 9 || *p == 13 || *p == 0 || *p == ',')) *token++ = *p++;

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
        if (*p == 0 || *p == 13) continue;

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
// loadDistribution() - Loads the fragment distribution definitions into RAM
//
// On Exit: the global "distribuitionList" object contains distribution definitions
//=================================================================================================
void loadDistribution()
{
    char fragmentName[1000];
    distribution_t distRecord;
    string line;
 
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
        if (*p == 0 || *p == 13) continue;

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

        // Clear the vector that will hold fragement names
        distRecord.fragments.clear();

        // Point to the comma separated fragement ids that come after the '$' delimeter
        p = delimeter;

        // Fetch every integer value after the name
        while (getNextCommaSeparatedToken(p, fragmentName))
        {
            if (fragment.find(fragmentName) == fragment.end())
            {
                throwRuntime("Undefined fragment name '%s'", fragmentName);
            }

            distRecord.fragments.push_back(fragmentName);
        }

        // And add this distribution record to the distribution list
        distributionList.push_back(distRecord);
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
        printf("%i : %i : %i   ", r.first, r.last, r.step);
        for (auto& f : r.fragments) printf(">%s<  ", f.c_str());
        printf("\n");

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

    // Loop through every record in the distribution list
    for (auto& distRec : distributionList)
    {
        // We're going to compute the length of this sequence of fragments
        uint32_t thisLength = 0;

        // Sum up the aggregate length of the fragments in this sequence
        for (auto& fragName : distRec.fragments)
        {
            thisLength += fragment[fragName].size();
        }

        // Keep track of the length of the longest sequence of fragments we find
        if (thisLength > longestLength) longestLength = thisLength;
    }

    // Hand the caller the length of the longest sequence of fragments
    return longestLength;
}
//=================================================================================================


void verifyDistributionIsValid()
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
}