/* 
 * filename: celltype.h
 * --------------------
 * This file defines a celltype class which will be 
 * instantiated as different cell types. At this point
 * it is not clear whether this should be a base class 
 * which will be inherited by unique cell classes or
 * whether all cells are instantiations of this class.  
 *
 */

#ifndef CELL_H_
#define CELL_H_

#include "cellparams.h"

class CellType
{
public:
	CellType(); 
	CellType(CellParams &cellParams);

	int getNumCellsX() const;
	int getNumCellsY() const;
	int getNumTotCells() const;

	std::string toString() const;

	// TODO: add an ID field so that we don't just compare numbers of cells...
	friend bool operator==(const CellType &thisCell, const CellType &otherCell);
	friend bool operator!=(const CellType &thisCell, const CellType &otherCell);
	friend std::ostream &operator<<(std::ostream &os, CellType &cell);

private:
	int numCellsX;
	int numCellsY;
	int numTotCells;
};

#endif /* CELL_H_ */
