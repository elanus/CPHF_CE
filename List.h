#ifndef ListH
#define ListH


#include <stdio.h>
#include "Element.h"
class List {
private:
public:
    Element* head;
    Element* tail;
    int size;
    List();
    ~List();
    void add(Element* ele);
    void remove(Element* ele);
    bool inList(int* p_data);
        bool inList(int a, int b, int c);
    string toString();
    Element* match1From(Element* start, int p1);
    int matchSofTCaller(int* sArray, int s, int t);
    int matchSofT(Element *ele, int* sArray, int s, int t, int i, int j);
};

#endif