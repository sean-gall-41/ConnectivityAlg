#ifndef DYNAMIC2DARRAY_
#define DYNAMIC2DARRAY_ 

#include <vector>
#include <cstddef>

template <typename type>
class Dynamic2DArray
{
public:
	Dynamic2DArray();
	Dynamic2DArray(int numRow, int numCol);
	Dynamic2DArray(int numRow, int numCol, type value);
	Dynamic2DArray(const Dynamic2DArray &copyFrom);
	~Dynamic2DArray();

	void fill(type value);
	int size() const;
	bool empty() const;
	void clear();	
	void reshape(int newNumRows, int newNumCols);
	type* data();

	// much easier thatn overloading '[]'
	type &operator()(int row, int col);	

	// equality	
	template <typename t>
	friend bool operator==(const Dynamic2DArray<t> &thisArr,
		const Dynamic2DArray<t> &otherArr);
	template <typename t>
	friend bool operator!=(const Dynamic2DArray<t> &thisArr,
		const Dynamic2DArray<t> &otherArr);

	Dynamic2DArray<type> &operator=(const Dynamic2DArray<type> &otherArr);

private:
	int numRows;
	int numCols;
	
	// in future, revert back to array
	std::vector<type> array;
};

#include "dynamic2darray_impl.h"

#endif /* DYNAMIC2DARRAY_ */ 

