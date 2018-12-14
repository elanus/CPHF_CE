//
//  Element.cpp
//  CPHF Optimizer & Verifier
//
//  Created by Erin Lanus on 1/29/15.
//  Copyright (c) 2015 Erin Lanus. All rights reserved.
//
#ifndef ElementPtrCPP
#define ElementPtrCPP

#include <stdio.h>
#include <iostream>
#include "Element.h"

using namespace std;

class ElementPtr {
private:
public:
    Element* ele;
    ElementPtr* prev;
    ElementPtr* next;
    
    ElementPtr(Element* p_ele)
    {
        ele = p_ele;
        prev = nullptr;
        next = nullptr;
    }
    
    ~ElementPtr()
    {
    }
    
    string toString()
    {
        return ele->toString();
    }
};

class ListPtr {
private:
public:
    ElementPtr* head;
    ElementPtr* tail;
    int size;
    ListPtr()
    {
        head = nullptr;
        tail = nullptr;
        size = 0;
    }
    
    ~ListPtr()
    {
        ElementPtr* curr = head;
        if(curr == nullptr)
            delete curr;
        else
        {
            while(curr->next != nullptr)
            {
                curr = curr->next;
                ElementPtr *temp = curr->prev;
                delete temp;
            }
            delete curr;
        }
    }
    void add(ElementPtr* ele)
    {
        if(head == nullptr)
        {
            head = ele;
            tail = ele;
        }
        else
        {
            tail->next = ele;
            ele->prev = tail;
            tail = ele;
        }
        size++;
    }
    string toString()
    {
        if(head == nullptr)
            return "()";
        ElementPtr* ele = head;
        string temp = "";
        while(ele != nullptr)
        {
            temp += ele->toString() + " ";
            ele = ele->next;
        }
        return temp;
    }
};

#endif