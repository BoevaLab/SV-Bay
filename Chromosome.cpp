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

#include "Chromosome.h"

using namespace std ;

Chromosome::Chromosome(int windowSize, std::string const& chromosome) {
	windowSize_=windowSize;
	chromosomeLengthBins_=0;
	someName_=chromosome;
}

Chromosome::~Chromosome(void)
{
	GenomeElements_.clear();
}

void Chromosome::addElement(std::string const& chromosome,float notNpercentage,int left, int right, int windowSize, std::string const& dnaSeq) {
	GenomeElements_.push_back(GenomeElement(chromosome, notNpercentage, left, right, windowSize, dnaSeq));
	chromosomeLengthBins_++;
}

void Chromosome::addElement(GenomeElement newGenomeElement) {
	GenomeElements_.push_back(newGenomeElement);
	chromosomeLengthBins_++;
}

int Chromosome::getLengthBins() {
	return chromosomeLengthBins_;
}

int Chromosome::findUntouched() {
	int guess = rand() % (chromosomeLengthBins_-2)+1;
	int count = 0;
	while (GenomeElements_[guess].isUntouchible()) {
		guess = rand() % (chromosomeLengthBins_-2)+1;
		count++;
		if (count == 10000) {
			guess=-1;
			break;
		}
	}
	return guess;
}

int Chromosome::find2Untouched() {
	if (chromosomeLengthBins_==0) {
		cout << "Chromosome "<< someName_ <<" of length 0!!!\n";
		cout << "Cannot find untouched element\n";
		return -1;
	}
	int guess = rand() % chromosomeLengthBins_;
	int count = 0;
	while (GenomeElements_[guess].isUntouchible()||GenomeElements_[guess+1].isUntouchible()) {
		guess = rand() % chromosomeLengthBins_;
		count++;
		if (count == 10000) {
			guess=-1;
			break;
		}
	}
	return guess;
}

void Chromosome::makeSmallInversion (int dupLength) {
	int selectedFragment = findUntouched();
	if (selectedFragment>=0)
		GenomeElements_[selectedFragment].makeSmallInversion(dupLength);
}

void Chromosome::makeLargeInversion(int numberOfElementsToInvert) {
	int selectedFragment = findUntouched();
	if (selectedFragment<0)
		return;
	int count = 0;
	while ( ((selectedFragment+numberOfElementsToInvert) > chromosomeLengthBins_) || GenomeElements_[selectedFragment+numberOfElementsToInvert-1].isUntouchible() ) {
		selectedFragment = findUntouched();
		if (selectedFragment<0)
			return;
		count ++;
		if (count >= 100000) {
			cout << "Impossible to find a place for a inversion of size  "<< numberOfElementsToInvert*windowSize_ << endl;
			return;
		}
	}

	vector<GenomeElement>::const_iterator first = GenomeElements_.begin() + selectedFragment;
	vector<GenomeElement>::const_iterator last = GenomeElements_.begin() + (selectedFragment+numberOfElementsToInvert);
	vector<GenomeElement> invertedRegion(first, last);
	std::reverse(invertedRegion.begin(),invertedRegion.end());
	for (std::vector<GenomeElement>::iterator it=invertedRegion.begin(); it!=invertedRegion.end(); ++it) {
		(*it).invert();
	}
	GenomeElements_.erase (GenomeElements_.begin()+selectedFragment,GenomeElements_.begin()+selectedFragment+ numberOfElementsToInvert);
	GenomeElements_.insert(GenomeElements_.begin()+selectedFragment,invertedRegion.begin(),invertedRegion.end());

	ostringstream convert;
	convert << "Inversion\t" << numberOfElementsToInvert*windowSize_ << "\t" << GenomeElements_[selectedFragment+numberOfElementsToInvert-1].getRightRefPos()
				<< "\t"<< GenomeElements_[selectedFragment].getLeftRefPos() << "\t" << GenomeElements_[selectedFragment+numberOfElementsToInvert-1].getRightRefSeq(FLANKSEQ)
				<< "\t"<< GenomeElements_[selectedFragment].getLeftRefSeq(FLANKSEQ) ;

	GenomeElements_[selectedFragment].setInsertedSV(convert.str());
	GenomeElements_[selectedFragment+1].setInsertedSV(convert.str());
}

void Chromosome::makeSmallDuplication (int dupLength) {
	int selectedFragment = findUntouched();
	if (selectedFragment>=0)
		GenomeElements_[selectedFragment].makeSmallDuplication(dupLength);
}

void Chromosome::makeSmallInsertion (int insLength) {
	int selectedFragment = findUntouched();
	if (selectedFragment>=0)
		GenomeElements_[selectedFragment].makeSmallInsertion(insLength);
}

void Chromosome::makeSmallDuplicationWithInversion (int dupLength) {
	int selectedFragment = findUntouched();
	if (selectedFragment>=0)
		GenomeElements_[selectedFragment].makeSmallDuplicationWithInversion(dupLength);
}
void Chromosome::makeLargeDuplicationWithInversion (int numberOfElementsToDuplicate){
	int selectedFragment = find2Untouched();
	if (selectedFragment<0)
			return;
	int count = 0;
	while ( ((selectedFragment+numberOfElementsToDuplicate) > chromosomeLengthBins_) || GenomeElements_[selectedFragment+numberOfElementsToDuplicate-1].isUntouchible() ) {
		selectedFragment = findUntouched();
		if (selectedFragment<0)
			return;
		count ++;
		if (count >= 100000) {
			cout << "Impossible to find a place for a duplication with inversion of size  "<< numberOfElementsToDuplicate*windowSize_ << endl;
			return;
		}
	}
	chromosomeLengthBins_ += numberOfElementsToDuplicate;

	vector<GenomeElement>::const_iterator first = GenomeElements_.begin() + selectedFragment;
	vector<GenomeElement>::const_iterator last = GenomeElements_.begin() + (selectedFragment+numberOfElementsToDuplicate);
	vector<GenomeElement> duplicatedRegion(first, last);

	ostringstream convert;
	convert << "Duplication_with_inversion\t" << numberOfElementsToDuplicate*windowSize_ << "\t" << duplicatedRegion[0].getLeftRefPos()
		<< "\t"<< duplicatedRegion.back().getRightRefPos() << "\t" << duplicatedRegion[0].getLeftRefSeq(FLANKSEQ)
				<< "\t"<< duplicatedRegion.back().getRightRefSeq(FLANKSEQ) ;

	duplicatedRegion[0].setInsertedSV(convert.str());
	duplicatedRegion.back().setInsertedSV(convert.str());

	std::reverse(duplicatedRegion.begin(),duplicatedRegion.end());
	for (std::vector<GenomeElement>::iterator it=duplicatedRegion.begin(); it!=duplicatedRegion.end(); ++it) {
			(*it).invert();
	}
	int coin = rand()%2;
	if (coin) {
		GenomeElements_.insert(GenomeElements_.begin()+selectedFragment,duplicatedRegion.begin(),duplicatedRegion.end()); //insert the inverted fragment before
	} else {
		GenomeElements_.insert(GenomeElements_.begin()+selectedFragment+numberOfElementsToDuplicate,duplicatedRegion.begin(),duplicatedRegion.end());//insert the inverted fragment after
	}
}


void Chromosome::makeFragmentInsertion(GenomeElement * NewElement, bool isInverted) {
	int selectedFragment = find2Untouched();
	ostringstream convert;
	int leftPosRef = GenomeElements_[selectedFragment].getRightRefPos() ;
	int rightPosRef = GenomeElements_[selectedFragment+1].getLeftRefPos();
	int fragmentLength = NewElement->getWindowSize();

	if(isInverted) {
		if(!NewElement->isInverted()) {
			NewElement->invert();
		}
	} else {
		if(NewElement->isInverted()) {
			NewElement->uninvert();
		}
	}

	string SVtype="Linking insertion\t";
	if (NewElement->isInverted()) {
		SVtype="Linking insertion with inversion\t";
	}

	convert << SVtype << fragmentLength << "\t" << leftPosRef
				<< "\t"<< rightPosRef << "\t" << GenomeElements_[selectedFragment].getRightRefSeq(FLANKSEQ)
				<< "\t"<< GenomeElements_[selectedFragment+1].getLeftRefSeq(FLANKSEQ) ;

	GenomeElements_[selectedFragment].setInsertedSV(convert.str());
	GenomeElements_[selectedFragment+1].setInsertedSV(convert.str());

	ostringstream convertIns;
	convertIns << SVtype << (*NewElement).getWindowSize() << "\t" << (*NewElement).getLeftRefPos()
				<< "\t"<< NewElement->getRightRefPos() << "\t" <<(*NewElement).getLeftRefSeq(FLANKSEQ)
				<< "\t"<< NewElement->getRightRefSeq(FLANKSEQ) ;
	NewElement->setInsertedSV(convertIns.str());
	GenomeElements_.insert((GenomeElements_.begin()+ selectedFragment + 1),(*NewElement));
	chromosomeLengthBins_++;
}

void Chromosome::amplifyMinute(int amplificationRate) {
	int unitLength=chromosomeLengthBins_;
	chromosomeLengthBins_ = unitLength*amplificationRate;

	vector<GenomeElement> duplicatedRegion(GenomeElements_.begin(), GenomeElements_.end());

	for (int i = 1; i < amplificationRate; i++) {
		//take the first unitLength elements and add them in the tail:
		GenomeElements_.insert(GenomeElements_.end(),duplicatedRegion.begin(),duplicatedRegion.end());
	}

	string type="Co-Amplification";
	if (unitLength==1) {type="Amplification";}

	for (int i = 0; i <unitLength; i++) {
		ostringstream convert;
		string isInvertedS="Direct";
		if (duplicatedRegion.at(i).isInverted()) {
			isInvertedS="Inverted";
		}
		convert << type<< " " << isInvertedS<<"\t" << duplicatedRegion.at(i).getWindowSize() << "\t" << duplicatedRegion.at(i).getLeftRefPos()
				<< "\t"<< duplicatedRegion.at(i).getRightRefPos() << "\t" << duplicatedRegion.at(i).getLeftRefSeq(FLANKSEQ)
				<< "\t"<< duplicatedRegion.at(i).getRightRefSeq(FLANKSEQ) ;
		GenomeElements_[i].setInsertedSV(convert.str());
	}
}

void Chromosome::makeLargeDuplication(int numberOfElementsToDuplicate) {
	int selectedFragment = findUntouched();
	if (selectedFragment<0)
			return;
	int count = 0;
	while ( ((selectedFragment+numberOfElementsToDuplicate) > chromosomeLengthBins_) || GenomeElements_[selectedFragment+numberOfElementsToDuplicate-1].isUntouchible() ) {
		selectedFragment = findUntouched();
		if (selectedFragment<0)
			return;
		count ++;
		if (count >= 100000) {
			cout << "Impossible to find a place for a duplication of size  "<< numberOfElementsToDuplicate*windowSize_ << endl;
			return;
		}
	}
	chromosomeLengthBins_ += numberOfElementsToDuplicate;

	vector<GenomeElement>::const_iterator first = GenomeElements_.begin() + selectedFragment;
	vector<GenomeElement>::const_iterator last = GenomeElements_.begin() + (selectedFragment+numberOfElementsToDuplicate);
	vector<GenomeElement> duplicatedRegion(first, last);

	GenomeElements_.insert(GenomeElements_.begin()+selectedFragment,duplicatedRegion.begin(),duplicatedRegion.end());

	ostringstream convert;
	convert << "Duplication\t" << numberOfElementsToDuplicate*windowSize_ << "\t" << GenomeElements_[selectedFragment].getLeftRefPos()
				<< "\t"<< GenomeElements_[selectedFragment+numberOfElementsToDuplicate-1].getRightRefPos() << "\t" << GenomeElements_[selectedFragment].getLeftRefSeq(FLANKSEQ)
				<< "\t"<< GenomeElements_[selectedFragment+numberOfElementsToDuplicate-1].getRightRefSeq(FLANKSEQ) ;

	GenomeElements_[selectedFragment].setInsertedSV(convert.str());
	GenomeElements_[selectedFragment+numberOfElementsToDuplicate-1].setInsertedSV(convert.str());
}

void Chromosome::makeSmallSizeDeletion (int delLength) {
	int selectedFragment = findUntouched();
	if (selectedFragment>=0)
		GenomeElements_[selectedFragment].makeSmallSizeDeletion(delLength);
}

void Chromosome::makeLargeDeletion(int numberOfElementsToDelete) {
	int selectedFragment = findUntouched();
	if (selectedFragment<0)
			return;
	int count = 0;
	while ( ((selectedFragment+numberOfElementsToDelete+2) > chromosomeLengthBins_) || GenomeElements_[selectedFragment+numberOfElementsToDelete+1].isUntouchible() ) {
		selectedFragment = findUntouched();
		if (selectedFragment<0)
			return;
		count ++;
		if (count >= 100000) {
			cout << "Impossible to find a place to delete  "<< numberOfElementsToDelete*windowSize_ << endl;
			return;
		}
	}
	chromosomeLengthBins_ -= numberOfElementsToDelete;
	GenomeElements_.erase (GenomeElements_.begin()+selectedFragment+1,GenomeElements_.begin()+selectedFragment+1+ numberOfElementsToDelete);

	ostringstream convert;
	convert << "Del\t" << numberOfElementsToDelete*windowSize_ << "\t" << GenomeElements_[selectedFragment].getRightRefPos()
				<< "\t"<< GenomeElements_[selectedFragment+1].getLeftRefPos() << "\t" << GenomeElements_[selectedFragment].getRightRefSeq(FLANKSEQ)
				<< "\t"<< GenomeElements_[selectedFragment].getLeftRefSeq(FLANKSEQ) ;

	GenomeElements_[selectedFragment].setInsertedSV(convert.str());
	GenomeElements_[selectedFragment+1].setInsertedSV(convert.str());

}

void Chromosome::makeTransclocation(Chromosome & chr1,Chromosome & chr2,bool isLeft1,bool isRight2) {
	int selectedFragment1 = chr1.find2Untouched();
	int selectedFragment2 = chr2.find2Untouched();
	if (selectedFragment1<0 || selectedFragment2<0)
		return;
	vector<GenomeElement> tail1,tail2;

	if(isLeft1) {
		vector<GenomeElement>::const_iterator first = chr1.GenomeElements_.begin() ;
		vector<GenomeElement>::const_iterator last = chr1.GenomeElements_.begin()+ selectedFragment1 + 1;
		tail1 = vector<GenomeElement> (first, last);
	} else {
		vector<GenomeElement>::const_iterator first = chr1.GenomeElements_.begin() + selectedFragment1 + 1;
		vector<GenomeElement>::const_iterator last = chr1.GenomeElements_.end();
		tail1 = vector<GenomeElement> (first, last);
		std::reverse(tail1.begin(),tail1.end());
		for (std::vector<GenomeElement>::iterator it=tail1.begin(); it!=tail1.end(); ++it) {
			(*it).invert();
		}
	}

	if(isRight2) {
		vector<GenomeElement>::const_iterator first = chr2.GenomeElements_.begin() + selectedFragment2 + 1;
		vector<GenomeElement>::const_iterator last = chr2.GenomeElements_.end();
		tail2 = vector<GenomeElement> (first, last);
	} else {
		vector<GenomeElement>::const_iterator first = chr2.GenomeElements_.begin() ;
		vector<GenomeElement>::const_iterator last = chr2.GenomeElements_.begin()+ selectedFragment2 + 1;
		tail2 = vector<GenomeElement> (first, last);
		std::reverse(tail2.begin(),tail2.end());
		for (std::vector<GenomeElement>::iterator it=tail2.begin(); it!=tail2.end(); ++it) {
			(*it).invert();
		}
	}

	int RefPos1, RefPos2;
	string seqFlank1, seqFlank2;

	RefPos1=tail1.back().getRightRefPos();
	seqFlank1 = tail1.back().getRightRefSeq(FLANKSEQ);
	string chr1Name = tail1.back().getChromosome();
	RefPos2=tail2[0].getLeftRefPos();
	seqFlank2 = tail2[0].getLeftRefSeq(FLANKSEQ);
	string chr2Name = tail2[0].getChromosome();

	ostringstream convert;
	convert << "Translocation\t" << chr1Name+chr2Name << "\t" << RefPos1
				<< "\t"<< RefPos2 << "\t" << seqFlank1
				<< "\t"<< seqFlank2 ;


	tail1.back().setInsertedSV(convert.str());
	tail2[0].setInsertedSV(convert.str());


	//connect the two tails:
	tail1.insert(tail1.end(), tail2.begin(), tail2.end());

	GenomeElements_ = tail1;
	chromosomeLengthBins_=GenomeElements_.size();

}

std::string Chromosome::getChromosomename() {
	return someName_;
}

void Chromosome::setChromosomename (std::string name) {
	someName_=name;
}

void Chromosome::printDonorFasta(std::ofstream & file,std::ofstream & commentsFile,int ID) {
	file << ">" << ID<<"_"<< someName_ << "\n";
	for (std::vector<GenomeElement>::iterator it=GenomeElements_.begin() ; it != GenomeElements_.end(); it++ ) {
		(*it).printDonorFasta(file,commentsFile,ID);
	}
}

GenomeElement * Chromosome::extractLargeFragment (int numberOfElementsToDuplicate) {
	if (chromosomeLengthBins_<numberOfElementsToDuplicate+10) {
		cout << "This chromosome "<<someName_<< " is too short to select a fragment of length " << numberOfElementsToDuplicate*windowSize_ <<endl;
		return NULL;
	}

	int selectedFragment = find2Untouched()+1;
	if (selectedFragment<0)
			return NULL;
	int count = 0;
	while ( ((selectedFragment+numberOfElementsToDuplicate) > chromosomeLengthBins_) || GenomeElements_[selectedFragment+numberOfElementsToDuplicate-1].isUntouchible() ) {
		selectedFragment = find2Untouched()+1;
		if (selectedFragment<0)
			return NULL;
		count ++;
		if (count >= 100000) {
			cout << "Impossible to find a fragment of size  "<< numberOfElementsToDuplicate*windowSize_ << endl;
			return NULL;
		}
	}

	vector<GenomeElement>::const_iterator first = GenomeElements_.begin() + selectedFragment;
	vector<GenomeElement>::const_iterator last = GenomeElements_.begin() + (selectedFragment+numberOfElementsToDuplicate);
	vector<GenomeElement> duplicatedRegion(first, last);

	if (duplicatedRegion.at(0).isInverted()) {  //TODO: should actually check isInverted for each segment
		std::reverse(duplicatedRegion.begin(),duplicatedRegion.end());
		for (std::vector<GenomeElement>::iterator it=duplicatedRegion.begin(); it!=duplicatedRegion.end(); ++it) {
				(*it).uninvert();
		}
	}

	float notNpercentage = 0;
	string dnaSeq="";
	for (int i=0; i<numberOfElementsToDuplicate; i++) {
		notNpercentage+=duplicatedRegion.at(i).getNotNpercentage();
		dnaSeq+=duplicatedRegion.at(i).getDnaSeq();
	}
	notNpercentage/=numberOfElementsToDuplicate;
	int left = duplicatedRegion.at(0).getLeftRefPos();
	int right=duplicatedRegion.at(0).getRightRefPos();
	GenomeElement * NewElement = new GenomeElement(duplicatedRegion.at(0).getChromosome(), notNpercentage, left, right, numberOfElementsToDuplicate*windowSize_, dnaSeq);

	ostringstream convert;
	convert << "Linking re-insertion (deleted part)\t" << numberOfElementsToDuplicate*windowSize_ << "\t" << left
				<< "\t"<< right << "\t" << NewElement->getLeftRefSeq(FLANKSEQ)
				<< "\t"<< NewElement->getRightRefSeq(FLANKSEQ) ;
	GenomeElements_.erase(GenomeElements_.begin()+selectedFragment,GenomeElements_.begin()+selectedFragment+numberOfElementsToDuplicate);
	GenomeElements_[selectedFragment-1].setInsertedSV(convert.str());
	GenomeElements_[selectedFragment].setInsertedSV(convert.str());
	chromosomeLengthBins_-=numberOfElementsToDuplicate;
	return NewElement;
}

GenomeElement * Chromosome::getLargeFragment (int numberOfElementsToDuplicate) {
	if (chromosomeLengthBins_<numberOfElementsToDuplicate+2) {
		cout << "This chromosome "<<someName_<< " is too short to select a fragment of length " << numberOfElementsToDuplicate*windowSize_ <<endl;
		return NULL;
	}

	int selectedFragment = findUntouched();
	if (selectedFragment<0)
			return NULL;
	int count = 0;
	while ( ((selectedFragment+numberOfElementsToDuplicate) > chromosomeLengthBins_) || GenomeElements_[selectedFragment+numberOfElementsToDuplicate-1].isUntouchible() ) {
		selectedFragment = findUntouched();
		if (selectedFragment<0)
			return NULL;
		count ++;
		if (count >= 100000) {
			cout << "Impossible to find a fragment of size  "<< numberOfElementsToDuplicate*windowSize_ << endl;
			return NULL;
		}
	}

	vector<GenomeElement>::const_iterator first = GenomeElements_.begin() + selectedFragment;
	vector<GenomeElement>::const_iterator last = GenomeElements_.begin() + (selectedFragment+numberOfElementsToDuplicate);
	vector<GenomeElement> duplicatedRegion(first, last);
	float notNpercentage = 0;
	string dnaSeq="";
	if (duplicatedRegion.at(0).isInverted()) {  //TODO: should actually check isInverted for each segment
		std::reverse(duplicatedRegion.begin(),duplicatedRegion.end());
		for (std::vector<GenomeElement>::iterator it=duplicatedRegion.begin(); it!=duplicatedRegion.end(); ++it) {
				(*it).uninvert();
		}
	}

	for (int i=0; i<numberOfElementsToDuplicate; i++) {
			notNpercentage+=duplicatedRegion.at(i).getNotNpercentage();
			dnaSeq+=duplicatedRegion.at(i).getDnaSeq();
	}


	notNpercentage/=numberOfElementsToDuplicate;
	GenomeElement * NewElement = new GenomeElement(duplicatedRegion.at(0).getChromosome(), notNpercentage, duplicatedRegion.at(0).getLeftRefPos(), duplicatedRegion.back().getRightRefPos(), numberOfElementsToDuplicate*windowSize_, dnaSeq);
	return NewElement;

}

GenomeElement * Chromosome::extractSmallFragment(int elementSize) {
	for (int i=0; i<15; i++) {
		int selectedFragment = find2Untouched();
		if (selectedFragment<0)
			return NULL;
		std::string const& chrName = GenomeElements_[selectedFragment].getChromosome();
		float notNpercentage = GenomeElements_[selectedFragment].getNotNpercentage();
		int ws = GenomeElements_[selectedFragment].getWindowSize();
		if (elementSize>ws) {
			cout << "Unpossible to extract a region of length "<< elementSize << " out of a fragment of length "<< ws << endl;
			cout << "Try N "<<i<<endl;
		}else {
			int startPos = rand() % (ws-elementSize);
			int left = GenomeElements_[selectedFragment].getLeftRefPos();
			int right = GenomeElements_[selectedFragment].getRightRefPos();
			if (left < right) {
				left = left+startPos;
				right = left+elementSize;
			}
			else {
				left = right+startPos;
				right = left+elementSize;
			}
			string dnaSeq = GenomeElements_[selectedFragment].getDnaSeq();
			dnaSeq=dnaSeq.substr(startPos,elementSize);
			GenomeElement * NewElement = new GenomeElement(chrName, notNpercentage, left, right, elementSize, dnaSeq);

//remove this fragment from the genome:
			dnaSeq = GenomeElements_[selectedFragment].getDnaSeq();
			dnaSeq = dnaSeq.substr(0,startPos)+dnaSeq.substr(startPos+elementSize);
			GenomeElements_[selectedFragment].setDNAsequence(dnaSeq);

			ostringstream convert;
			convert << "Linking re-insertion (deleted region)\t" << elementSize << "\t" << left
				<< "\t"<< right << "\t" << NewElement->getLeftRefSeq(FLANKSEQ)
				<< "\t"<< NewElement->getRightRefSeq(FLANKSEQ) ;

			GenomeElements_[selectedFragment].setInsertedSV(convert.str());
			return NewElement;
		}
	}
	return NULL;
}


GenomeElement * Chromosome::getSmallFragment (int elementSize) {

	for (int i=0; i<15; i++) {
		int selectedFragment = findUntouched();
		if (selectedFragment<0)
			return NULL;
		std::string const& chrName = GenomeElements_[selectedFragment].getChromosome();
		float notNpercentage = GenomeElements_[selectedFragment].getNotNpercentage();
		int ws = GenomeElements_[selectedFragment].getWindowSize();
		if (elementSize>ws) {
			cout << "Unpossible to extract a region of length "<< elementSize << " out of a fragment of length "<< ws << endl;
			cout << "Try N "<<i<<endl;
		}else {
			int startPos = rand() % (ws-elementSize);
			int left = GenomeElements_[selectedFragment].getLeftRefPos();
			int right = GenomeElements_[selectedFragment].getRightRefPos();
			if (left < right) {
				left = left+startPos;
				right = left+elementSize;
			}
			else {
				left = right+startPos;
				right = left+elementSize;
			}
			string dnaSeq = GenomeElements_[selectedFragment].getDnaSeq();
			dnaSeq=dnaSeq.substr(startPos,elementSize);
			GenomeElement * NewElement = new GenomeElement(chrName, notNpercentage, left, right, elementSize, dnaSeq);
			return NewElement;
		}
	}
	return NULL;
}
