#pragma once
#include <stddef.h>

// Source : https://stackoverflow.com/questions/12693263/std-deque-is-surprisingly-slow

template <class T>
class FastStack {

    public:
        T* st;
        int allocationSize;
        int lastIndex;

        FastStack(int stackSize);
        FastStack();
        ~FastStack();

        inline void resize(int newSize);
        inline void push(T x);
        inline void pop();
        inline T getAndRemove();
        inline T getLast();
        inline void clear();
};

template <class T>
FastStack<T>::FastStack() {

    lastIndex = -1;
    st = NULL;
}

template <class T>
FastStack<T>::FastStack(int stackSize) {

    st = NULL;
    this->allocationSize = stackSize;
    st = new T[stackSize];
    lastIndex = -1;
}

template <class T>
FastStack<T>::~FastStack() {

    delete [] st;
}

template <class T>
void FastStack<T>::clear() {

    lastIndex = -1;
}

template <class T>
T FastStack<T>::getLast() {

    return st[lastIndex];
}

template <class T>
T FastStack<T>::getAndRemove() {

    return st[lastIndex--];
}

template <class T>
void FastStack<T>::pop() {

    --lastIndex;
}

template <class T>
void FastStack<T>::push(T x) {

    st[++lastIndex] = x;
}

template <class T>
void FastStack<T>::resize(int newSize) {

    if (st != NULL) delete [] st;
    st = new T[newSize];
}
