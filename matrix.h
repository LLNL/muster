#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
using namespace std;

/**
 * Simple class for matrices with static dimensions.
 * TODO: add bounds checking.
 */
template <class T> 
class matrix {
private:

    static const int DEFAULT_ROWS = 2;
    static const int DEFAULT_COLS = 2;

    /** Data in the matrix. */
    vector< vector<T> > data;
    typedef typename vector< vector<T> >::iterator vIterator;

public:

    /** 
     * Constructs a new dynamic matrix, optionally reserving space enough for
     * the dimensions specified.
     * 
     * @param rows number of rows in matrix
     * @param cols number of cols in matrix
     */
    matrix(const unsigned _rows = DEFAULT_ROWS, const unsigned _cols = DEFAULT_COLS) :
        data(_rows, vector<T>(_cols))
    { }


    /** Free up memory for this matrix. */
    virtual ~matrix() {
        //nothing necessary.
    }

    /** 
     * Reserves space for specified number of rows and cols, if there's not that
     * much already.
     */
    void reserve(unsigned rows, unsigned cols) {
        data.reserve(rows);
        for (vIterator v=data.begin(); v != data.end(); v++) {
            (*v).reserve(cols);
        }
    }


    /** retrieves a reference to the element (col, row) in the matrix. */
    T& get(const unsigned row, const unsigned col) {
        return data[row][col];
    }

    /** retrieves a reference to the element (col, row) in the matrix. */
    const T& get(const unsigned row, const unsigned col) const {
        return data[row][col];
    }

    /** Shorthand for get() */
    T& operator()(const unsigned row, const unsigned col) {
        return data[row][col];
    }

    /** Shorthand for get() */
    const T& operator()(const unsigned row, const unsigned col) const {
        return get(row, col);
    }

    /** 
     * Takes two vectors of indices, and returns a new 
     * rowindices.size() by colindices.size() matrix where the (i,j)th 
     * elements in this matrix are at row rowindices[i] and column colindices[j]
     * in this matrix.
     */
    template <class InputIterator>
    matrix<T> *getSubmatrix(
        InputIterator rowStart, InputIterator rowEnd,
        InputIterator colStart, InputIterator colEnd
    ) const {
        matrix<T> *m = new matrix<T>(0,0);

        vector<T> v;  //new vector to copy into each row
        int row = 0;

        for (InputIterator r = rowStart; r != rowEnd; r++) {
            m->data.push_back(v);
            vector<T>& rv = m->data[row];  //ref to copied vector in new matrix
            for (InputIterator c = colStart; c != colEnd; c++) {
                rv.push_back(get(*r, *c));
            }
            row++;
        }
        return m;
    }

        
    /** Prints contents out, if possible. */
    void printOn(ostream& out = cout) const {
        for (unsigned r=0; r < data.size(); r++) {
            for (unsigned c = 0; c < data[r].size(); c++) {
                out << data[r][c] << " ";
            }
            out << endl;
        }
    }

};


#endif //MATRIX_H
