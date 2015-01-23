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
#ifndef _CHR_GENOMEEL_H
#define _CHR_GENOMEEL_H

#include <vector>
#include <map>
#include <stdlib.h>
#include <iostream>

#include "myFunc.h"


class GenomeElement
{
public:
    GenomeElement(std::string const& chrName);
	GenomeElement(std::string const& chrName, float notNpercentage, int left, int right, int windowSize, std::string const& dnaSeq);
	GenomeElement(std::string const& chrName, float notNpercentage, int left, int right, int windowSize, std::string const& dnaSeq,std::vector <int> delPos,std::vector <int> delSizes, std::vector <int> gainPos,std::vector <int> gainSizes);
	~GenomeElement(void);

	std::string getChromosome();
	int getRightRefPos() ;
	int getLeftRefPos();
	int getWindowSize();
	std::string getRightRefSeq(int length);
	std::string getLeftRefSeq(int length);
	std::string getDnaSeq();
	float getNotNpercentage();

	bool isUntouchible();

	void makeSmallSizeDeletion(int delLength);
	void makeSmallDuplication(int dupLength);
	void makeSmallInversion(int dupLength);
	void makeSmallInsertion(int insLength);
	void makeSmallDuplicationWithInversion(int dupLength);
	void setInsertedSV (std::string const& message);
	void setDNAsequence(std::string const& dnaSeq);

	void printDonorFasta(std::ofstream & file,std::ofstream & commentsFile,int ID);
	void invert();
	void uninvert();
	bool isInverted();

private:
    int leftRefPos_;
	int rightRefPos_;
	bool isInverted_;
	std::string chromosome_;
	std::string dnaSeq_;
	/*std::vector <int> gains_;
	std::vector <int> dels_;	*/
	std::string InsertedSVs_;
	float notNpercentage_;
	int windowSize_;
};
#endif

