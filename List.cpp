//
//  List.cpp
//  CPHF Optimizer & Verifier
//
//  Created by Erin Lanus on 1/29/15.
//  Copyright (c) 2015 Erin Lanus. All rights reserved.
//
#ifndef Listcpp
#define Listcpp

#include "List.h"


List::List()
{
    head = nullptr;
    tail = nullptr;
    size = 0;
}
List::~List()
{
    Element* curr = head;
    if(curr == nullptr)
        delete curr;
    else
    {
        while(curr->next != nullptr)
        {
            curr = curr->next;
            Element *temp = curr->prev;
            //curr->prev->~Element();
            //temp->~Element();
            delete temp;
        }
        //curr->~Element();
        delete curr;
    }
    
}
void List::add(Element* ele)
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
void List::remove(Element* ele)
{
    if(head == ele && tail == ele)
    {
        head = nullptr;
        tail = nullptr;
    }
    else if(head==ele)
    {
        head = ele->next;
        head->prev = nullptr;
    }
    else if(tail == ele)
    {
        tail = ele->prev ;
        tail->next = nullptr;
    }
    else
    {
        ele->prev->next = ele->next;
        ele->next->prev = ele->prev;
    }
    delete ele;
    size--;
}
string List::toString()
{
    if(size ==0)
        return "()";
    Element* ele = head;
    string temp = "";
    while(ele != nullptr)
    {
        temp += ele->toString() + " ";
        ele = ele->next;
    }
    return temp;
}

bool List::inList(int a, int b, int c)
{
    cout << "using wrong version of inList";
    return false;
}

bool List::inList(int* p_data)
{
    bool found = false;
    Element* curr = head;
    while(curr != nullptr && !found)
    {
        bool same = true;
        for(int i = 0; i < Element::t && same; i++)
        {
            if(curr->data[i] != p_data[i])
                same = false;
        }
        if(same)
            found = true;
        curr = curr->next;
    }
    return found;
}

// Examine the list for elements containing one column starting at a certain element in the list until found
Element* List::match1From(Element* start, int p1)
{
    Element *curr = start;
    while(curr != nullptr)
    {
        bool match = false;
        for(int i = 0; i < Element::t && !match; i++)
        {
            if(curr->data[i] == p1)
                match = true;
        }
        if(match)
            break;
        else
            curr = curr->next;
    }
    return curr;
}

int List::matchSofT(Element *ele, int* sArray, int s, int t, int i, int j)
{
    int val = 0;
    if(j == s)
    {
        return 1;
    }
    if(ele->data[i] == sArray[j])
    {
        val = matchSofT(ele, sArray, s, t, i+1, j+1);
    }
    else if(i < t + j - s)
    {
        val = matchSofT(ele, sArray, s, t, i+1, j);
    }
    return val;
}
#endif