#ifndef ElementH
#define ElementH

#include <stdio.h>
#include <iostream>

using namespace std;

class Element {
private:
public:
    int *data;
    static int t;
    Element* prev;
    Element* next;
    
    Element(int *p_data);
    //~Element();
    string toString();
};
#endif