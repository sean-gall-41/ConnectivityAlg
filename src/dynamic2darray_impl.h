#include <algorithm>

const int DEFAULT_SIZE = 1;

template<typename type>
Dynamic2DArray<type>::Dynamic2DArray()
{
	numRows     = DEFAULT_SIZE;
	numCols     = DEFAULT_SIZE;
	
	array.resize(numRows * numCols);
}

template <typename type>
Dynamic2DArray<type>::Dynamic2DArray(int numRows, int numCols):
	numRows(numRows), numCols(numCols)
{
	array.resize(numRows * numCols);	
}

template <typename type>
Dynamic2DArray<type>::Dynamic2DArray(int numRows, int numCols, type value):
	numRows(numRows), numCols(numCols)
{
	array.resize(numRows * numCols);
	fill(value);
}

template <typename type>
Dynamic2DArray<type>::Dynamic2DArray(const Dynamic2DArray &copyFrom)
{
	Dynamic2DArray<type>(copyFrom.numRows, copyFrom.numCols);	
	std::copy(copyFrom.array.begin(), copyFrom.array.end(), this->array.begin());
}

template <typename type>
Dynamic2DArray<type>::~Dynamic2DArray() {}

template <typename type>
void Dynamic2DArray<type>::fill(type value)
{
	std::fill(array.begin(), array.end(), value);	
}

template <typename type>
int Dynamic2DArray<type>::size() const
{
	return array.size();
}

template <typename type>
bool Dynamic2DArray<type>::empty() const
{
	return array.empty();
}

template <typename type>
void Dynamic2DArray<type>::clear() 
{
	array.clear();
}

template <typename type>
void Dynamic2DArray<type>::reshape(int newNumRows, int newNumCols)
{

	array.resize(newNumRows * newNumCols);
	
	numRows = newNumRows;
	numCols = newNumCols;

}

template <typename type>
type &Dynamic2DArray<type>::operator()(int row, int col)
{
	return array[row * numCols + col];
}	

template <typename t>
bool operator==(const Dynamic2DArray<t> &thisArr,
	const Dynamic2DArray<t> &otherArr)
{ 
	return (thisArr.numRows == otherArr.numRows)
		&& (thisArr.numCols == otherArr.numCols) 
		&& (thisArr.array == otherArr.array);
}

template <typename t>
bool operator!=(const Dynamic2DArray<t> &thisArr,
	const Dynamic2DArray<t> &otherArr) 
{
	return !(thisArr == otherArr);
}

template <typename type>
Dynamic2DArray<type> &Dynamic2DArray<type>::operator=(const Dynamic2DArray<type> &otherArr)
{
	if (this != &otherArr)
	{
		this->array   = otherArr.array;
		this->numRows = otherArr.numRows;
		this->numCols = otherArr.numCols;
	}

	return *this;
}

