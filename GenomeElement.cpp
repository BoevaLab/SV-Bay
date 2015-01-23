/*************************************************************************
Copyright (c) 2014-2016, Valentina BOEVA.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/


#include "GenomeElement.h"

using namespace std ;

GenomeElement::GenomeElement(std::string const& chrName) {
	chromosome_ = chrName;
	leftRefPos_=0;
	rightRefPos_=0;
	isInverted_=0;
	windowSize_=0;
	std::string dnaSeq_="";
	InsertedSVs_="";
	notNpercentage_=0;
}

GenomeElement::GenomeElement(std::string const& chrName, float notNpercentage, int left, int right, int windowSize, std::string const& dnaSeq) {
	chromosome_ = chrName;
	leftRefPos_=left;
	rightRefPos_=right;
	isInverted_=0;
	windowSize_=windowSize;
	dnaSeq_=dnaSeq;
	//gains_.resize(windowSize,0);
	//dels_.resize(windowSize,0);
	InsertedSVs_="";
	notNpercentage_=notNpercentage;
}

GenomeElement::GenomeElement(std::string const& chrName, float notNpercentage, int left, int right, int windowSize, std::string const& dnaSeq, std::vector <int> delPos,std::vector <int> delSizes, std::vector <int> gainPos,std::vector <int> gainSizes) {
	chromosome_ = chrName;
	leftRefPos_=left;
	rightRefPos_=right;
	windowSize_=windowSize;
	dnaSeq_=dnaSeq;
	//gains_.resize(windowSize,0);
	//dels_.resize(windowSize,0);
	notNpercentage_=notNpercentage;

	//for(int i=0;i<delPos.size();i++) {
	//	dels_[delPos[i]-leftRefPos_]=delSizes[i];
	//}
	//for(int i=0;i<gainPos.size();i++) {
	//	dels_[gainPos[i]-leftRefPos_]=gainSizes[i];
	//}
	InsertedSVs_="";
}

std::string GenomeElement::getChromosome() {
	return chromosome_;
}

GenomeElement::~GenomeElement(void)
{
	//dels_.clear();
	//gains_.clear();
	dnaSeq_ = "";
}

int GenomeElement::getWindowSize() {
 return windowSize_;
}

bool GenomeElement::isUntouchible() {
	if (InsertedSVs_.length()>0 || notNpercentage_ <= 0.5) {
		return 1;
	}
	else return 0;
}

float GenomeElement::getNotNpercentage() {
	return notNpercentage_;
}


void GenomeElement::makeSmallInsertion(int insLength) {
	int startingPos = rand() % (windowSize_ - FLANKSEQ*2-1)+FLANKSEQ;
	ostringstream convert;
	convert << "Insertion\t" << insLength<< "\t" << (leftRefPos_+startingPos)<< "\t"<< (leftRefPos_+startingPos+insLength)<< "\t"
		<< dnaSeq_.substr(startingPos-FLANKSEQ,FLANKSEQ) << "\t"<< dnaSeq_.substr(startingPos,FLANKSEQ);
	InsertedSVs_= convert.str();
	string insertedRegion = randomSeq(insLength);
	dnaSeq_.insert(startingPos, insertedRegion);
}

void GenomeElement::makeSmallInversion(int invLength){
	int startingPos = rand() % (windowSize_ - invLength-1);
	string initialRegion = dnaSeq_.substr(startingPos,invLength);
	string invertedRegion = complementary(initialRegion);
	dnaSeq_.erase(startingPos, invLength);
	dnaSeq_.insert(startingPos, invertedRegion);
	windowSize_=dnaSeq_.length();
	int leftFlankLength=min(FLANKSEQ,invLength);
	int rightFlankLength=leftFlankLength;
	ostringstream convert;
	convert << "Inversion\t" << invLength<< "\t" << (leftRefPos_+startingPos)<< "\t"<< (leftRefPos_+startingPos+invLength)<< "\t"
		<< initialRegion.substr(0,leftFlankLength) << "\t"<< initialRegion.substr(invLength-rightFlankLength);

	InsertedSVs_= convert.str();
}

void GenomeElement::makeSmallDuplication(int dupLength) {
	int startingPos = rand() % (windowSize_ - dupLength-1);
	string duplicatedRegion = dnaSeq_.substr(startingPos,dupLength);
	dnaSeq_.insert(startingPos, duplicatedRegion);
	windowSize_=dnaSeq_.length();
	int leftFlankLength=min(FLANKSEQ,dupLength);
	int rightFlankLength=leftFlankLength;
	ostringstream convert;
	convert << "Duplication\t" << dupLength<< "\t" << (leftRefPos_+startingPos)<< "\t"<< (leftRefPos_+startingPos+dupLength)<< "\t" << duplicatedRegion.substr(0,leftFlankLength) << "\t"<< duplicatedRegion.substr(dupLength-rightFlankLength);
	InsertedSVs_= convert.str();
}

void GenomeElement::makeSmallDuplicationWithInversion(int dupLength) {
	if (windowSize_ - dupLength -1 <=0 ) {
		cout << "..Unable to duplicate "<< dupLength <<"bp\n";
		return;
	}
	int startingPos = rand() % (windowSize_ - dupLength -1)+1;
	string duplicatedRegion = dnaSeq_.substr(startingPos,dupLength);
	string duplicatedRegionInverted = complementary(duplicatedRegion);

	int coin = rand()%2;
	if (coin) {
		dnaSeq_.insert(startingPos, duplicatedRegionInverted);
	} else {
		dnaSeq_.insert(startingPos+dupLength, duplicatedRegionInverted);
	}

	windowSize_=dnaSeq_.length();
	int leftFlankLength=min(FLANKSEQ,dupLength);
	int rightFlankLength=leftFlankLength;
	ostringstream convert;
	convert << "Duplication\t" << dupLength<< "\t" << (leftRefPos_+startingPos)<< "\t"<< (leftRefPos_+startingPos+dupLength)<< "\t" << duplicatedRegion.substr(0,leftFlankLength) << "\t"<< duplicatedRegion.substr(dupLength-rightFlankLength);
	InsertedSVs_= convert.str();
}

void GenomeElement::makeSmallSizeDeletion(int delLength) {
	int startingPos = rand() % (windowSize_ - delLength-1);

	dnaSeq_.erase(startingPos, delLength);
	windowSize_=dnaSeq_.length();
	int leftFlankLength=min(FLANKSEQ,startingPos);
	int rightFlankLength=min(FLANKSEQ,windowSize_-startingPos-1);
	ostringstream convert;
	convert << "Del\t" << delLength<< "\t" << (leftRefPos_+startingPos)<< "\t"<< (leftRefPos_+startingPos+delLength)<< "\t" << dnaSeq_.substr(startingPos-leftFlankLength,leftFlankLength) << "\t"<< dnaSeq_.substr(startingPos,rightFlankLength);
	InsertedSVs_= convert.str();
}

void GenomeElement::setInsertedSV (std::string const& message) {
	InsertedSVs_=message;
}

void GenomeElement::setDNAsequence(std::string const& dnaSeq) {
	dnaSeq_=dnaSeq;
	windowSize_=dnaSeq.length();
}

int GenomeElement::getRightRefPos() {
	if (!isInverted_) return rightRefPos_;
	return leftRefPos_;
}
int GenomeElement::getLeftRefPos() {
	if (!isInverted_) return leftRefPos_;
	return rightRefPos_;
}
std::string GenomeElement::getDnaSeq() {
	return dnaSeq_;
}

std::string GenomeElement::getLeftRefSeq(int length) {
	if (!isInverted_) return dnaSeq_.substr(windowSize_-FLANKSEQ,FLANKSEQ);
	return dnaSeq_.substr(0,FLANKSEQ);
}
std::string GenomeElement::getRightRefSeq(int length) {
	if (!isInverted_) return dnaSeq_.substr(0,FLANKSEQ);
	return dnaSeq_.substr(windowSize_-FLANKSEQ,FLANKSEQ);
}

void GenomeElement::printDonorFasta(std::ofstream & file, std::ofstream & commentsFile,int ID) {
	string toPrint;
	if (!isInverted_)
		toPrint = dnaSeq_;
	else
		toPrint = complementary(dnaSeq_);
	for (int i = 0; i < windowSize_; i+=FASTASTRINGLEN) {
		file << toPrint.substr(i,FASTASTRINGLEN) << "\n";
	}
	if (InsertedSVs_.length()>0)
		commentsFile << ID<< "_" << chromosome_<<"\t"<< InsertedSVs_ << "\n";
}

void GenomeElement::invert() {
	isInverted_ = 1;
}

void GenomeElement::uninvert() {
	isInverted_ = 0;
}

bool GenomeElement::isInverted () {
	return isInverted_;
}
