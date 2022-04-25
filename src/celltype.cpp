#include <cassert>
#include "celltype.h"

CellType::CellType() : numCellsX(0), numCellsY(0), numTotCells(0) {};

CellType::CellType(CellParams &cellParams)
{
	numCellsX = std::stoi(cellParams["numCellsX"]);
	numCellsY = std::stoi(cellParams["numCellsY"]);
	numTotCells = std::stoi(cellParams["numTotCells"]);
	assert(numTotCells == numCellsX * numCellsY);
}

int CellType::getNumCellsX() const {return numCellsX;}

int CellType::getNumCellsY() const {return numCellsY;}

int CellType::getNumTotCells() const {return numTotCells;}

std::string CellType::toString() const
{
	std::string outString = "";
	outString += "[ numCellsX : " + std::to_string(numCellsX) + " ]\n";
	outString += "[ numCellsY : " + std::to_string(numCellsY) + " ]\n";
	outString += "[ numTotCells: " + std::to_string(numTotCells) + " ]";
	return outString;
}

bool operator==(const CellType &thisCell, const CellType &otherCell)
{
	return (thisCell.numCellsX == otherCell.numCellsX)
		&& (thisCell.numCellsY == otherCell.numCellsY)
		&& (thisCell.numTotCells == otherCell.numTotCells);	
}

bool operator!=(const CellType &thisCell, const CellType &otherCell)
{
	return !(thisCell == otherCell);
}

std::ostream &operator<<(std::ostream &os, CellType &cell)
{
	return os << cell.toString();
}

