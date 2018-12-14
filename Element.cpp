//
//  Element.cpp
//  CPHF Optimizer & Verifier
//
//  Created by Erin Lanus on 1/29/15.
//  Copyright (c) 2015 Erin Lanus. All rights reserved.
//
#ifndef ElementCPP
#define ElementCPP

#include "Element.h"

int Element::t = 3;

Element::Element(int *p_data)
{
    data = p_data;
    prev = nullptr;
    next = nullptr;
}

//Element::~Element()
//{
//    delete data;
//}

string Element::toString()
{
    string temp = "(";
    for(int i = 0; i < t-1; i++)
    {
        temp += to_string(data[i]) + ", ";
    }
    temp += to_string(data[t-1]) + ")";
    return temp;
}

#endif