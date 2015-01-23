/*************************************************************************
Copyright (c) 2014-2016, Valentina BOEVA.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#pragma once
#ifndef _CHR_LIST_H
#define _CHR_LIST_H

#include "GenomeElement.h"

class Chromosome
{
public:
    Chromosome(int windowSize,std::string const& chromosome);
	~Chromosome(void);
	void addElement(std::string const& chromosome,float notNpercentage,int left, int right, int windowSize, std::string const& text);
	void addElement(GenomeElement newGenomeElement);

	int getLengthBins();
	std::string getChromosomename();
	void setChromosomename (std::string name);

	GenomeElement * getSmallFragment (int elementSize); //will keep this fragment in the initial place in the genome;
	GenomeElement * getLargeFragment (int numberOfElementsToDuplicate); //will keep this fragment in the initial place in the genome;

	GenomeElement * extractSmallFragment (int elementSize); //will remove this fragment from the initial place in the genome
	GenomeElement * extractLargeFragment (int numberOfElementsToDuplicate); //will remove this fragment from the initial place in the genome

	int findUntouched();
	int find2Untouched();

	void makeSmallSizeDeletion (int delLength) ;
	void makeLargeDeletion(int numberOfElementsToDelete);
	void makeSmallDuplication (int dupLength) ;
	void makeLargeDuplication(int numberOfElementsToDuplicate);
	void makeSmallDuplicationWithInversion (int dupLength) ;
	void makeLargeDuplicationWithInversion (int numberOfElementsToDuplicate);
	void makeSmallInversion (int dupLength) ;
	void makeLargeInversion(int numberOfElementsToDuplicate);
	void makeTransclocation(Chromosome & chr1,Chromosome & chr2,bool isLeft1,bool isRight2);
	void makeSmallInsertion (int insLength) ;
	void makeFragmentInsertion(GenomeElement * NewElement, bool isInverted);
	void amplifyMinute(int amplificationRate);

	void printDonorFasta(std::ofstream & file,std::ofstream & commentsFile,int ID);

private:

    std::vector<GenomeElement> GenomeElements_;
	int windowSize_;
	int chromosomeLengthBins_;
	std::string someName_;
};
#endif
