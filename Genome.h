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
#ifndef _GENOME_CPN_H
#define _GENOME_CPN_H

#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <fstream>

#include "tgsim.h"
#include "Chromosome.h"
#include "myFunc.h"
//#include <cstdint> //for windows
#include <stdint.h>  //for Linux

class Genome
{
public:
	Genome(int windowSize);
	Genome(int windowSize, std::string const& genomeFastaFilename);
	~Genome(void);

	int getNumberOfChromosomes();

	void addFromFasta(std::string const& genomeFastaFilename);

	int getWindowSize(void);
	void addDeletions(int number, int minLength, int maxLength);
	void addDuplications(int number, int minLength, int maxLength);
	void addDuplicationsWithInversion(int number, int minLength, int maxLength);
	void addInversions(int number, int minLength, int maxLength);
	void addTranslocations(int number, bool isBalanced);
	void addUnbalancedTranslocations(int number);
	void addSmallInsertions(int number, int minLength, int maxLength);
	void addLinkingInsertions(int number, int minLength, int maxLength, bool isInverted);
	void addReInsertions(int number, int minLength, int maxLength, bool isInverted);
	void addCoAmplification(int number, int minLength, int maxLength,bool sameChr,int amplificationRate);

	void deleteChromosomes(int number);
	void duplicateChromosomes(int number);
	void duplicateGenome();

	void printDonorGenome(std::string const& fastaFileName, std::string const& commentsFileName);
	int findChromForSV(int SVlength);

private:
	int windowSize_;
	std::vector<Chromosome> chromosomes_;
};
#endif

