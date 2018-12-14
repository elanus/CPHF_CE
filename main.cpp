/*  Author: Erin Lanus
    Affiliation: Arizona State University
 
    Implements a Conditional Expectation Algorithm to build Covering Perfect Hash Families
    Version works for parameters:
                t = 3, q <= 7, q a prime power
                t = 4, q <= 16, q a prime power
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <time.h>
#include "List.h"
#include "ElementPtr.cpp"

using namespace std;

bool debug, verbose;
ofstream output;
string basePath = "./";
int k, q, t, N, qtminus1;
int row; // the current row of the SCPHF being built
int** combinations;
bool* ncMatrix;
long nct = 0; // number of non-covering tuples, dependent upon t and q)
double P = 0; // probably of a t-subset being covered if all columns free and selected randomly, for conditional expectation/density
double threshold = P;
int earlyExitCount = 0;
int maxRows = 10; // number of rows in storage for SCPHF, start with 10, double when needed
time_t start, finish; // used for clocking time to run
int currNonCov = 0; // number of elements in the noncovering list so far

int primepowList[5] = {4, 8, 9, 16, 25};
int** addTable;
int** mulTable;
int** subTable;
int** subPrimeTable; // used for constructing subTable for prime powers

long nctTable3[] = {0,0,0,12,80,300,0,1960};
long nctTable4[] = {0,0,0,117936,3636864,46593000,0,2012501736, 8905430016, 32702473440, 0, 298290835920, 0, 1872470514096, 0, 0, 18323862036480}; // upperbound for q > 7
long *nctTable[5];

// creates a table of combinations (i choose j) up to (k choose t)
void buildCombiTable()
{
    combinations = new int*[k+1];
    // base cases are 0 and 1
    for(int i = 0; i<2; i++)
    {
        combinations[i] = new int[t+1];
        for(int j=0; j<t+1; j++)
            combinations[i][j] = (j == 0 || j == i);
    }
    // rest are defined by recursion
    for(int i = 2; i<k+1; i++)
    {
        combinations[i] = new int[t+1];
        combinations[i][0] = 1;
        for(int j = 1; j<t+1; j++)
        {
            combinations[i][j] = combinations[i-1][j] + combinations[i-1][j-1];
        }
    }
    // print
    if(debug)
    {
        cout << "Printing the combinations table: " << endl;
        for(int i=0;i<k+1;i++)
        {
            for(int j=0;j<t+1; j++)
                cout << combinations[i][j] << "\t";
            cout << endl;
        }
        cout << endl;
    }
}

// checks to see if the value given is a prime power up to the value allowed by the program
bool primepowCheck(int n)
{
    bool flag = false;
    for(int i = 0; i < 5 && !flag; i++)
    {
        if(primepowList[i] == n)
            flag = true;
    }
    return flag;
}

// Creates the arithmetic tables for operating in the field
// If q is a prime, we compute the tables directly
// If q is one of the first 5 prime powers, we read the tables from a file
void setupArithmeticTables()
{
    // checks to see if q is one of the first 5 prime powers (up to and including 25)
    // if yes, we read in the addition, mult, sub tables for v
    bool isprimepower = primepowCheck(q);
    if(isprimepower)
    {
        if(debug)
            cout << "v is one of the first 5 powers (exp > 1) of a prime " << endl;
        
        ifstream pptables;
        string pptablespath = basePath + "/primepowertables.txt";
        pptables.open(pptablespath);
        
        if(pptables.is_open())
        {
            string op;
            int val, curr_q;
            bool found_q = false;
            while(!found_q)
            {
                pptables >> curr_q;
                cout << "curr_q " << curr_q << endl;
                // if that's the table we need to read, then...
                if(curr_q == q)
                {
                    found_q = true;
                    pptables >> op;
                    addTable = new int*[q];
                    for(int i = 0; i < q; i++) // for each element of GF(q) read in,
                    {
                        addTable[i] = new int[q]; // we make a row in the addition table that holds q elements
                        for(int j = 0; j < q; j++) // and read in the elements for that row
                        {
                            pptables >> val;
                            addTable[i][j] = val;
                        }
                    }
                    pptables >> op;
                    mulTable = new int*[q];
                    for(int i = 0; i < q; i++) // for each element of GF(q) read in,
                    {
                        mulTable[i] = new int[q]; // we make a row in the multiplicatino table that holds q elements
                        for(int j = 0; j < q; j++) // and read in the elements for that row
                        {
                            pptables >> val;
                            mulTable[i][j] = val;
                        }
                    }
                    pptables >> op;
                    subPrimeTable = new int*[q];
                    for(int i = 0; i < q; i++) // for each element of GF(q) read in,
                    {
                        subPrimeTable[i] = new int[q]; // we make a row in the subPrime table that holds q elements
                        for(int j = 0; j < q; j++) // and read in the elements for that row
                        {
                            pptables >> val;
                            subPrimeTable[i][j] = val;
                        }
                    }
                    
                }
                // if it's not what we need to read, skip through so many lines...
                else
                {
                    for(int j = 0; j < 3; j++)  // 3 sets of tables, add, mul, sub to skip over
                    {
                        pptables >> op;
                        for(int i = 0; i < curr_q*curr_q; i++)
                            pptables >> val;
                    }
                }
            }
            pptables.close();
        }
        else
        {
            cout << "Error opening primepowertables. Exiting..." << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    // If not a prime power, we assume it's prime, so compute tables using mod
    // We do not handle incorrect values for v
    else
    {
        // +
        addTable = new int*[q];
        for(int i = 0; i < q; i++)
        {
            addTable[i] = new int[q];
            for(int j = 0; j < q; j++)
                addTable[i][j] = (i + j) % q;
        }
        // *
        mulTable = new int*[q];
        for(int i = 0; i < q; i++)
        {
            mulTable[i] = new int[q];
            for(int j = 0; j < q; j++)
                mulTable[i][j] = (i * j) % q;
        }
    }
    // -
    // used for component-wise subtraction for permutation vector headers, that is 020 - 020 = 000; 011 - 020 = 021, etc
    // computes x - y in the field
    // e.g. if x = x_0q^0 + x_1q^1 + ... + x_(t-2)q^(t-2) and y = y_0q^0 + y_1q^1 + ... + y_(t-2)q^(t-2)
    // only difference in prime vs prime power is that prime does coordinate level subtraction directly while prime power looks it up from the table we imported
    // the subtraction table is an q^(t-1) square array so we have two indices, x and y, with the same range
    subTable = new int*[qtminus1];
    for(int x = 0; x < qtminus1; x++)
    {
        subTable[x] = new int[qtminus1];
        for(int y = 0; y < qtminus1; y++)
        {
            // compute the values for the vectors: each position is the remainder of the quotient / q and the new quotient carries through to the next position
            int quotientx = x;
            int quotienty = y;
            int remainderx, remaindery;
            int value = 0;
            for(int pos = 0; pos < t-1; pos++)
            {
                // compute xbaseq
                remainderx = quotientx % q;
                quotientx = quotientx / q;
                //compute ybaseq
                remaindery = quotienty % q;
                quotienty = quotienty / q;
                
                if(isprimepower)
                    value += subPrimeTable[remainderx][remaindery] * pow(q,pos);
                else
                {
                    if(remainderx - remaindery < 0) // can't operate on negative numbers so roll it around mod q
                        remainderx += q;
                    value += (remainderx - remaindery) * pow(q,pos);
                }
            }
            subTable[x][y] = value;
        }
    }
    // Print the tables
    if(debug)
    {
        cout << "+" << endl;
        for(int i = 0; i < q; i++)
        {
            for(int j = 0; j < q; j++)
                cout << addTable[i][j] << "\t";
            cout << endl;
        }
        
        cout << "*" << endl;
        for(int i = 0; i < q; i++)
        {
            for(int j = 0; j < q; j++)
                cout << mulTable[i][j] << "\t";
            cout << endl;
        }
        
        cout << "-" << endl;
        for(int i = 0; i < qtminus1; i++)
        {
            for(int j = 0; j < qtminus1; j++)
                cout << subTable[i][j] << "\t";
            cout << endl;
        }
    }
    // add the counts of noncovering tuples to the table
    for(int i = 0; i < 3; i++)
        nctTable[i] = new long[0];
    nctTable[3] = nctTable3;
    nctTable[4] = nctTable4;
}

// copy the values in tuple1 into tuple2
void copyTuple(int* tuple1, int* tuple2, int length)
{
    for(int i = 0; i < length; i++)
        tuple2[i] = tuple1[i];
}

// print a tuple to std out
void printTuple(int* tuple, int length)
{
    for(int i = 0; i < length; i++)
        cout << tuple[i] << "\t";
}
// sorts the symbols of a t-tuple into increasing order in place in the int* tuple parameter
void sortTuple(int* tuple, int length)
{
    for(int a = 1; a < length; a++)
    {
        int b = a;
        int c = b-1;
        while(c >= 0)
        {
            if(tuple[b] < tuple[c])
            {
                int temp = tuple[c];
                tuple[c] = tuple[b];
                tuple[b] = temp;
                b -= 1;
                c -= 1;
            }
            else
                break;
        }
    }
}

// returns the tuple of columns in set array of length size that correspond to rank
int* rankToTuple(int* set, int rank, int size)
{
    if(size==0)
        return set;
    
    int m = 0;
    while(combinations[m+1][size] <= rank)
        m+=1;
    int newRank = rank - combinations[m][size];
    set = rankToTuple(set,newRank,size-1);
    set[size-1]= m; // set ct = m
    return set;
}

// Creates a new SCPHF with twice as many rows, copies the values, and frees the old SCPHF
int** doubleArray(int** Array, int rows, int cols)
{
    int** temp = new int*[rows*2];
    for(int i = 0; i < rows; i++)
    {
        temp[i] = new int[cols];
        for(int j = 0; j < cols; j++)
            temp[i][j] = Array[i][j];
    }
    for(int i = rows; i < rows*2; i++)
    {
        temp[i] = new int[cols];
        for(int j = 0; j < cols; j++)
            temp[i][j] = -1;
    }
    // free memory for old array
    for(int i = 0; i < rows; i++)
        delete[] Array[i];
    delete[] Array;
    
    maxRows = rows*2;
    return temp;
}

// prints the SCPHF array
void printArray(int** array, int rows, int cols)
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols-1; j++)
        {
            cout << array[i][j] << "\t";
        }
        cout << array[i][cols-1] << endl;
    }
}

// given a sorted t-tuple, produce the integer index for the boolean ncMatrix
// computes the canonical form of the tuple by shifting e.g. if tuple is 1 2 6 8, produces tuple 1-1, 2-1, 6-1, 8-1 = 0 1 5 7
// then for the tuple 0 1 5 7, we would ignore the 0 and compute by taking the sum of (value of this symbol) * q^(t-1) ^ (number of positions from the end)
int computeNcMatrixIndex(int* tuple)
{
    int index = 0;
    // compute the reduced form by subtracting the first tuple element from all others (could have chosen any element to subtract from all)
    // and then finally subtract the first element from itself to make it 0
    for(int i = t-1; i >= 0; i--)
        tuple[i] = subTable[tuple[i]][tuple[0]];
    
    // sort the reduced form as they may no longer be in order to put in canonical form
    sortTuple(tuple, t);
    
    // compute the index of the canonical form
    for(int i = 1; i < t; i++)
    {
        // looks up the shifted value tuple[i] - tuple[0] with arithmetic in the field
        // and then multiplies by it's position to build the index one position at a time
        index += tuple[i] * pow(qtminus1, t-1-i);
    }
    
    // print the canonical form of the tuple and the index into the ncMatrix
    /*if(debug)
     {
     cout << "tuple in canonical form: ";
     for(int i = 0; i < t; i++)
        cout << tuple[i] << " ";
    cout << "\tindex " << index << endl;
     }*/
    
    return index;
}

// given a tuple, returns yes if the tuple is covering; returns no if the tuple is noncovering
// calls computeNcMatrixIndex on the tuple to get the index and then returns the opposite of the t/f value stored in the ncMatrix (matrix indicating noncovering tuples)
bool isCovering(int* tuple)
{
    return !ncMatrix[computeNcMatrixIndex(tuple)];
}
// returns true if all symbols in the tuple are distinct, returns false if any symbol appears more than once
bool isDistinct(int* tuple, int size)
{
    // must have no repeated symbol to have a hope of being covered
    bool allDistinct = true;
    for(int i = 1; i < size && allDistinct; i++)
    {
        if(tuple[i] == tuple[i-1])
            allDistinct = false;
    }
    return allDistinct;
}

// ncMatrix is a lookup to see if a tuple is noncovering.
// We convert tuples to a canonical form in which the first symbol is 0 so we only need t-1 indices to lookup whether the tuple is noncovering
// The memory allocation is t-dependent, so since we can't create an (t-1)-dimensional array at runtime, we create a 1-D array with q^(t-1)^(t-1) rows
// Then we read in the values from a file
bool* createNcMatrix()
{
    // the ncMatrix has dimension t-1 and each dimension has length q^(t-1)
    int matrixSize = pow(qtminus1,t-1);
    ncMatrix = (bool *) malloc (matrixSize * sizeof(bool));
    // initialize the matrix to all 0s, add 1s only for items read from file
    for(int i = 0; i < matrixSize; i++)
        ncMatrix[i] = 0;
    
    ifstream input;
    string NCpath = basePath + "NCReduced/noncovering_reduced_"+to_string(q)+"_"+to_string(t)+".txt";
    input.open(NCpath);
    if(input.is_open())
    {
        // for debugging, print the noncovering tuples and their indices into the ncMatrix
        if(debug)
            cout << "creating noncovering matrix" << endl;
        
        // count the number of distinct (unordered) noncovering tuples
        int distinctNoncov = 0;
        
        // create a t-1 lenth array to store the t values on each line
        int* data = new int[t];
        
        while(!input.eof())
        {
            // read in t values at a time. Assumes input file is correctly formatted
            // can code without data[] by reading into a variable and then using immediately, but this looks easier to follow
            for(int i = 0; i < t; i++)
                input >> data[i];
            
            // compute the index in ncMatrix of the noncovering tuple and set that position to true
            int index = computeNcMatrixIndex(data);
            ncMatrix[index] = 1;
            
            // increase the count of distinct noncovering tuples seen
            distinctNoncov++;
        }
        
        // for debugging, print the ncMatrix
        /*if(debug)
         {
         for(int i = 0; i < matrixSize; i++)
         cout << ncMatrix[i] << "\t";
         }*/
        
        cout << "count of noncovering list in reduced memory/canonical form: " << distinctNoncov << endl;
        nct = nctTable[t][q];
        input.close();
        delete[] data;
    }
    else
    {
        cout << "Could not read file for non-covering list located at " << NCpath << endl;
        cout << "Exiting..." << endl;
        exit(EXIT_FAILURE);
    }
    
    return ncMatrix;
}

// determines which t-subsets remain to be covered, control depends on 3 cases:
// if it's nonserial first row, then everything is remaining initially
// treat serial first row just like an addon array where there is just 1 row added (addedRows = 1)
// if adding on to a already created SCPHF, need to check all of the rows for each of the t-subsets, and add to remaining only if not covered by any row
List* createRemainingList(int** SCPHF, List* remaining, int level, int lowerbound)
{
    // use ranks to iterate through all (k choose t) t-subsets of columns in the SCPHF
    int kct = combinations[k][t];
    int* columns = (int*) malloc(t * sizeof(int));
    int* symbols = (int*) malloc(t * sizeof(int));
    for(int rank = 0; rank < kct; rank++)
    {
        // unrank to get the columns of the t-subset in colexicographical order
        columns = rankToTuple(columns, rank, t);
        
        // start by assuming these columns aren't covered
        bool covered = false;
        
        // if serial, addedRows = 1, so this does exactly the same thing
        // this is the opportunity to check if it has been covered; it's it's neither addOn or serial first row, we get past here with covered still = false
        if(row > 0)
        {
            for(int r = 0; r < row && !covered; r++)
            {
                // otherwise, gather the symbols by putting the symbol from the rth row of the SCPHF in the ith column into the ith symbol slot
                for(int i = 0; i < t; i++)
                    symbols[i]=SCPHF[r][columns[i]];
                
                // sort them
                sortTuple(symbols, t);
                
                //if first thing is -1, something not placed - everything in tuple has to be placed to be possibly covering
                if(symbols[0] != -1)
                {
                    // if the symbols are all distinct and it is a covering tuple, then mark it as covered
                    if(isDistinct(symbols,t) && isCovering(symbols))
                        covered = true;
                }
            }
        }
        // if it is not covered by some row of the array we are adding on to, add it to remaining list
        if(!covered)
        {
            int *data = (int*) malloc(t * sizeof(int));
            copyTuple(columns, data, t);
            remaining->add(new Element(data));
        }
    }
    free(columns);
    free(symbols);
    return remaining;
}

// recurisive function that returns the number of ways to complete a t-tuple with some symbols already fixed
// kill branches that result in tuples of non-distinct symbols
int computeCountCovering(int* partialTuple, int numFixed)
{
    // count will be an expectation
    double count = 0;
    // base case, all symbols are fixed, return 1 if the (now fully fixed) tuple is covering and 0 otherwise
    if(numFixed == t)
        count = isCovering(partialTuple);
    else
    {
        // use a copy of the partialTuple passed in as we will be sorting it in place and different calls to this function use the same stored partialTuple
        int* tempTuple = new int[t];
        for(int i = 0; i < t; i++)
            tempTuple[i] = -1;

        // look at every way to place a symbol in the first free cell
        // if the partial tuple is distinct, count the ways of finishing the tuple with that symbol that are covering
        for(int symbol = 0; symbol < qtminus1; symbol++)
        {
            copyTuple(partialTuple, tempTuple, t);
            // store the symbol in the numFixed position since array indices are off by one
            tempTuple[numFixed] = symbol;
            /*if(debug)
            {
                cout << "printing the temporary tuple ";
                printTuple(tempTuple, numFixed+1);
                cout << endl;
            }*/
            sortTuple(tempTuple, numFixed+1);
            if(isDistinct(tempTuple, numFixed+1))
                count += computeCountCovering(tempTuple, numFixed+1);
        }
        delete[] tempTuple;
    }
    return count;
}

// returns the symbol with the highest expectation for the cell at row row and column focusCol
int mostCoveringSymbol(List *remaining, int** SCPHF, int row, int focusCol, int defaultSym)
{
    int maxSymbol = defaultSym;
    if(debug)
        cout << "focusCol " << focusCol << endl;
    
    // copies tsubs from remainingtsubs into tsubsInvolvingCol to look ahead at interactions of focusCol with all tsubs
    ListPtr* tsubsInvolvingCol = new ListPtr();
    ListPtr* only1Free = new ListPtr();
    Element* currPos = remaining->head;
    // collect all uncovered tsubs involving focusCol that have at least 1 col fixed
    while(currPos != nullptr)
    {
        Element* found = remaining->match1From(currPos, focusCol);
        if(found != nullptr)
        {
            if(debug)
                cout << "found match " << found->toString() << endl;
            // check that at least one other col has a symbol fixed. If all are free, there are no restrictions and no symbols are better than others
            int numFixed = 0;
            for(int i = 0; i < t; i++)
            {
                if(SCPHF[row][found->data[i]] != -1)
                    numFixed++;
            }
            
            if(numFixed > 0)
                tsubsInvolvingCol->add(new ElementPtr(found));
            if((numFixed+1) == t)
                only1Free->add(new ElementPtr(found));
            currPos = found->next; // whether added or skipped, look for the next one
        }
        else
            currPos = found;
    }
    
    if(debug)
        cout << "tsubsInvolvingCol for " << focusCol << ": " << tsubsInvolvingCol->toString();
    
    int size = tsubsInvolvingCol->size;
    if(debug)
        cout << "size=" << size << endl;
    
    //if no remaining tsubs involve focusCol, just place the default symbol
    if(size == 0)
        SCPHF[row][focusCol] = defaultSym;
    else
    {
        // Considering each symbol, check all tsubs in which this col is involved and look for the symbol that covers the most
        double max = 0;
        // boolean flags used for early exit -- if a symbol is found that covers "enough", set early exit,
        bool earlyExit = false;
        int symbolsChecked = 0;
        if(debug)
            cout << "expected @ random: " << P * size << endl;
        // look through all of the symbols unless the early exit flag is set by reaching the threshhold
        for(int i = defaultSym; symbolsChecked < qtminus1 && !earlyExit; i = (i+1) % qtminus1)
        {
            symbolsChecked++;
            double expectation = 0;
            // look at each of the remaining t-subsets that involve the focus column
            ElementPtr* currPos = tsubsInvolvingCol->head;
            while(currPos != nullptr)
            {
                // have a tsub (currPos) - get the other symbols
                // numFixed is the number of fixed symbols - start with t and reduce it for every -1 seen
                int numFixed = t;
                
                // currpos->data is cols, let symbols be the symbols in those cols
                int* symbols = new int[t];

                for(int j = 0; j < t; j++)
                {
                    symbols[j] = SCPHF[row][currPos->ele->data[j]];
                    if(symbols[j] == -1)
                        numFixed--;
                }
                
                // collect the fixed symbols + i in an array and sort it
                int* sorted = new int[numFixed+1];
                int j = 0;
                for(int l = 0; l < t; l++)
                {
                    if(symbols[l] != -1)
                    {
                        sorted[j] = symbols[l];
                        j++;
                    }
                }
                sorted[numFixed] = i;
                sortTuple(sorted, numFixed+1);
                
                // if the current value being considered, i, isn't already placed in one of the other columns being considered,
                // or if the two other columns being considered don't already have the same value, continue
                if(isDistinct(sorted,numFixed+1))
                {
                    // computeCountCovering returns the number of ways sorted could be completed as a covering tuple
                    int count = computeCountCovering(sorted, numFixed+1);
                    // the expected number is the ways to complete it as a covering tuple / ways to complete it
                    // there are numFixed+1 symbols placed in sorted, so there are t-(numFixed+1) additional symbols to choose from the q^(t-1) symbols
                    expectation += (double) count / pow(qtminus1, t - (numFixed+1));
                }
                currPos = currPos->next;
                delete[] sorted;
                delete[] symbols;
            }
            if(debug)
                cout << " i = " << i << ", expectation = " << expectation << endl;
            
            if(expectation > max)
            {
                max = expectation;
                maxSymbol = i;
            }
            // early exit if something found that does at least as well as threshold % of size
            if(expectation >= threshold * size)
            {
                earlyExit = true;
                earlyExitCount++;
            }
        }
        
        SCPHF[row][focusCol] = maxSymbol;
        if(debug)
        {
            cout << "Printing t-subs with only1Free ";
            cout << only1Free->toString() << endl;
            
        }
        ElementPtr* curr = only1Free->head;
        while(curr != nullptr)
        {
            int* symbols = new int[t];
            for(int j = 0; j < t; j++)
                symbols[j] = SCPHF[row][curr->ele->data[j]];
            
            sortTuple(symbols, t);
            
            if(isDistinct(symbols,t) && isCovering(symbols))
            {
                if(debug)
                    cout << "removing " << curr->ele->toString() << " from remaining" << endl;
                remaining->remove(curr->ele); // delete the tsub of COLUMNS from remainingtsubs
            }
            delete[] symbols;
            
            curr = curr->next;
        }
    }
    only1Free->~ListPtr();
    delete tsubsInvolvingCol;
    return maxSymbol;
}

// reads a SCPHF from an input file
int** readFile()
{
    int** SCPHF;
    int fileRows = 0;
    
    ifstream input;
    string inputpath = basePath + "input.txt";
    input.open(inputpath);
    string line;
    
    // scans the input file for the number of columns and rows
    // update this section with an appropriately formatted file header so we can skip this section
    if(input.is_open())
    {
        cout << "Input should be space delimited. First line of input file should be: N k q t" << endl;
        input >> fileRows >> k >> q >> t;
        maxRows = fileRows;
        // creates a SCPHF in storage with enough columns and rows
        SCPHF = new int*[fileRows];
        // initialize cells of the array to the "unplaced" symbol, -1
        for(int r = 0; r < fileRows; r++)
        {
            SCPHF[r] = new int[k];
            for(int c = 0; c < k; c++)
                input >> SCPHF[r][c];
        }
    }
    input.close();
    
    // might as well double check the format
    input.open(inputpath);
    if(input.is_open())
    {
        string line;
        int countrow, countcol;
        
        getline(input, line);
        countrow = 1;
        
        while(getline(input, line))
        {
            countrow++;
            countcol = 0;
            for(int i=0; i<line.length(); i++)
                if(line[i] == ' ')
                    countcol++;
            if(countcol+1 != k)
            {
                cout << "Incorrect number of columns Exiting...";
                exit(EXIT_FAILURE);
            }
        }
        if(countrow != (fileRows+1))
        {
            cout << "Incorrect number of rows. Exiting...";
            exit(EXIT_FAILURE);
        }
    }
    input.close();
    
    // print the array
    for(int i = 0; i < fileRows; i++)
    {
        for(int j = 0; j < k; j++)
            cout << SCPHF[i][j] << " ";
        cout << endl;
    }
    row = fileRows;
    return SCPHF;
}

// builds a SCPHF one row at a time
int** buildSCPHF(int** SCPHF, bool serialFirst, bool addOn)
{
    List* remaining = new List();
    
    //cout << "Printing remaining tsubs before first row..." << endl;
    //cout << remaining->toString();
    int sym = 0;
    
    // serialFirst = use a completely deterministic 1st row - it's "as good" as any other option since more repeated symbols is always worse
    // and we could permute any set of columns using as few repeated symbols as possible to appear in "serial" form
    if(serialFirst)
    {
        for(int i = 0; i < k; i++)
            SCPHF[0][i] = i % qtminus1;
        
        row = 1;
    }
    
    createRemainingList(SCPHF, remaining, 0, 0);
    
    /* Row 2+ Construction --- make it a while loop to check if remainingtsubs is empty */
    
    /***************************** Make Additional Rows *****************************/
    /* master while loop that continues making rows until everything covered */
    while(remaining->size != 0)
    {
        if(row >= maxRows)
            SCPHF = doubleArray(SCPHF, maxRows, k);
        if(row == 0)
        {
            cout << "Initial number of tsubs " << remaining->size << endl;
            output << "Initial number of tsubs " << remaining->size << endl;
        }
        else
        {
            cout << "Number of tsubs left after row " << row-1 << ": " << remaining->size <<  endl;
            output << "Number of tsubs left after row " << row-1 << ": " << remaining->size <<  endl;
            if(debug)
                cout << remaining->toString() << endl;
        }

        // placing 0 and 1 in the first positions of a row does not change the density of the row
        SCPHF[row][0] = 0;
        SCPHF[row][1] = 1;
        
        //sym is the default or "starting" value, always reset sym to 2 at start of additional row
        sym = 2;
        
        // Density of a row is the number of expected t-subsets of columns to be covered for the first time by the row
        // e.g. the ceiling of the remaining t-subsets times the probability of choosing a covering tuple at random
        // then the expected remaining size is the current remaining size minus the density
        int density = ceil(remaining->size * P);
        int expecRemainingSize = (remaining->size - density);
        cout << "Density of row " << row << ": " << density << " , expected remaining after: " << expecRemainingSize << endl;
        output << "Density of row " << row << ": " << density << " , expected remaining after: " << expecRemainingSize << endl;
        
        // look at each free column in the array in the array
        for(int col = 2; col < k; col++)
        {
            if(SCPHF[row][col] == -1)
            {
                int symbol = mostCoveringSymbol(remaining, SCPHF, row, col, sym);
                // update sym to try to get as balanced set of sym in a row as possible
                if(symbol == sym)
                    sym = (sym + 1) % qtminus1; // move to next symbol in set
                if(debug)
                    cout << "symbol is " << symbol << " column is " << col << endl;
            }
        }
        
        if(debug)
        {
            cout<<endl;
            for(int i = 0; i < k; i++)
                cout << SCPHF[row][i] << "\t";
            cout << "\n\nRow " << row << " completed\n-----------------------------------------" << endl;
        }
        
        row++;
        
        if(remaining->size > expecRemainingSize)
            cerr << "MORE REMAINING THAN EXPECTED!!! " << remaining->size << " " << expecRemainingSize << endl;
        
    } // end of while loop that continues making rows until everything covered
    
    N = row;
    remaining->~List();
    return SCPHF;
}

// Takes in parameters for the SCPHF and reports final result as well as time to completion
int main(int argc, const char * argv[])
{
    int** SCPHF;
    bool serialFirst = false;
    bool addOn = false;
    string noyes = " (0=no, 1=yes) ";
    cout << "SCPHF Constructor by Conditional Expectation. Version for t > 2, with memory reduction implementation" << endl;
    cout << "-----------------------------------------------------------------------------------------------------" << endl;
    cout << "Print output for debugging?" << noyes;
    cin >> debug;
    cout << "Verbose mode? " << noyes;
    cin >> verbose;
    cout << "Add on to existing (from input)?" << noyes;
    cin >> addOn;
    if(addOn)
        SCPHF = readFile();
    else
    {
        cout << "Serial first row?";
        if(verbose)
            cout << "(useful when k large relative to number of symbols)";
        cout << noyes;
        cin >> serialFirst;
        
        cout << "Enter k: ";
        if(verbose)
            cout << "(number of columns) ";
        cin >> k;
        
        if(verbose)
        {
            cout << "This program works when q is prime or q is a prime power up to and including 25. If a non-valid value for q is entered, the program will not work as expected. Symbols in the SCPHF will be from set q^(t-1). ";
        }
        cout << "Enter q: ";
        cin >> q;
        
        cout << "Enter the strength, t > 2: ";
        cin >> t;
        if(t < 3)
        {
            cout << "Invalid value for t. Exiting..." << endl;
            exit(EXIT_FAILURE);
        }
        SCPHF = new int*[maxRows];
        // initialize cells of the array to the "unplaced" symbol, -1
        for(int r = 0; r < maxRows; r++)
        {
            SCPHF[r] = new int[k];
            for(int c = 0; c < k; c++)
                SCPHF[r][c] = -1;
        }
    }
    Element::t = t;
    
    time(&start);
    
    qtminus1 = (int)pow(q, t-1);
    setupArithmeticTables();
    buildCombiTable();
    createNcMatrix();
    
    long distinct = 1;
    for(int i = 0; i < t; i++)
        distinct *= (qtminus1-i);
    
    P = (distinct - nct)/pow(qtminus1,t);
    cout   << "number of ordered distinct t-tuples: " << distinct << "\nnumber of (ordered) noncovering tuples (nct): " << nct << endl;
    cout << "Number of ordered covering tuples, distinct - nct: " << (distinct - nct)  << "\nnumber of all tuples, (q^(t-1))^t " << pow(qtminus1,t) <<  endl;
    cout << "Probability of choosing a covering tuple at random, P: " << P << endl;
    cout << "\nEnter the threshold as percent of total to meet in search (enter value between 0 and 1) ";
    if(verbose)
        cout << "(choose threshold at least as large as P to meet bounds) ";
    cin >> threshold;
    cout << "threshold " << threshold << endl;
    time(&finish);
    double seconds = difftime(finish, start);
    int minutes = seconds / 60;
    seconds = seconds - (minutes * 60);
    printf("Setup time: %d minutes %.f seconds\n", minutes, seconds);
    
    string outputpath = basePath + "Results/";
    outputpath += "SCPHF(" + to_string(k) + "," + to_string(qtminus1) + "," + to_string(t) + ")";
    outputpath += to_string(threshold) + "_" + to_string(serialFirst) + "_" + to_string(addOn) + ".txt";
    output.open(outputpath, ofstream::app);
    if(output.is_open())
    {
        time(&start);
        output << "distinct t-tuples " << distinct << ", nct " << nct << ", distinct - nct " << (distinct - nct)  << ", pow(qtminus1,t) " << pow(qtminus1,t) <<  endl << "P " << P << endl;
        
        cout << "--------------------------------------------------------------------------------------------\n";
        output << "--------------------------------------------------------------------------------------------\n";
        cout   << "k=" << k << ", q=" << q << ", q^(t-1)=" << qtminus1 << ", t=" << t << ", serial first row " << serialFirst << ", threshold " << threshold << endl;
        output << "k=" << k << ", q=" << q << ", q^(t-1)=" << qtminus1 << ", t=" << t << ", serial first row " << serialFirst << ", threshold " << threshold << endl;
        cout << "--------------------------------------------------------------------------------------------\n";
        output << "--------------------------------------------------------------------------------------------\n";
        
        // passes the partial SCPHF either blank or read from file and returns the finished SCPHF
        SCPHF = buildSCPHF(SCPHF, serialFirst, addOn);
        
        cout   << "SCPHF(" << N << ";" << to_string(k) << "," << to_string(qtminus1) << "," << to_string(t) << ")" << endl;
        output << "SCPHF(" << N << ";" << to_string(k) << "," << to_string(qtminus1) << "," << to_string(t) << ")" << endl;
        for(int r = 0; r < N; r++)
        {
            for(int c = 0; c < k-1; c++)
            {
                cout   << SCPHF[r][c] << " ";
                output << SCPHF[r][c] << " ";
            }
            cout   << SCPHF[r][k-1] << endl;
            output << SCPHF[r][k-1] << endl;
        }
        cout   << "and produces a CA(" << (N * (pow(q, t) - q) + q) << ";" << t << "," << k << "," << q << ")" << endl;
        output << "and produces a CA(" << (N * (pow(q, t) - q) + q) << ";" << t << "," << k << "," << q << ")" << endl;
        
        cout << "Early exits " << earlyExitCount << endl;
        time(&finish);
        double seconds = difftime(finish, start);
        int minutes = seconds / 60;
        seconds = seconds - (minutes * 60);
        cout << minutes << " minutes " << seconds << " seconds" << endl;
        output << minutes << " minutes " << seconds << " seconds" << endl;
        output.close();
    }
    else
    {
        cout << "Problem opening output file.";
        exit(EXIT_FAILURE);
    }
    
    // clean up
    for(int i = 0; i < maxRows; i++)
        delete[] SCPHF[i];
    delete[] SCPHF;
    for(int i = 0; i < k+1; i++)
        delete[] combinations[i];
    free(ncMatrix);
    for(int i = 0; i < q; i++)
    {
        delete[] addTable[i];
        delete[] mulTable[i];
    }
    delete[] addTable;
    delete[] mulTable;
    return 0;
}
