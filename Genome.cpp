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


#include "Genome.h"

using namespace std ;

Genome::Genome(int windowSize)
{
	windowSize_=windowSize;
}

Genome::Genome(int windowSize, std::string const& genomeFastaFilename) {
	cout << "Reading the genome from "<< genomeFastaFilename <<endl;
	windowSize_=windowSize;
	ifstream file;
	file.open(genomeFastaFilename.c_str());
	if(!file.is_open()) {
		file.clear();
		cerr << genomeFastaFilename <<"cannot be opened\n";
		exit (-1);
	}

	int chromosomeIndex = 0;
	int startOnChr = 0;
	char letter; // here we read the chromosome letter by letter. Can do it better!
	int count = 0;
	int countN = 0;
	string line;
    string text = "";
	string chromosome = "";

	while (!file.eof()) {
		file>>letter;
		if (letter == '>') {
			if (chromosome!="") { //save the previous segment if you want
				float notNpercentage = float(count-countN)/count;
				int left = startOnChr;
				int lastFragmentLength=text.length();
				startOnChr+=lastFragmentLength;
				chromosomes_.back().addElement(chromosome,notNpercentage,left,startOnChr-1,lastFragmentLength,text);
			}
			getline (file,chromosome);
			cout << "reading chromosome "<< chromosome << "\n";
			Chromosome newChromosome(windowSize,chromosome);
			chromosomes_.push_back(newChromosome);
			//reset:
			countN = 0;
			count = 0;
			text="";
			startOnChr = 0;
		} else {
			countN += isN(letter);
			count ++;
			text+=letter;  //text.push_back(letter);
			if (count == windowSize_) {
				float notNpercentage = float(count-countN)/count;
				int left = startOnChr;
				startOnChr+=windowSize_;
				chromosomes_.back().addElement(chromosome,notNpercentage,left,startOnChr-1,windowSize_,text);
				//reset:
				countN = 0;
				count = 0;
				text="";
			}
		}
	}
	file.close();
	cout << "End-of-file reached.." << endl;
	/*vector <int> smallChromosomes;
	for (int i=0; i< chromosomes_.size();i++) {
		if (chromosomes_.at(i).getLengthBins()==0) {
			smallChromosomes.push_back(i);
		}
	}
	std::reverse(smallChromosomes.begin(),smallChromosomes.end());
	for (int i=0; i<smallChromosomes.size(); i++) {
		cout << "..Removing chr " << chromosomes_.at(smallChromosomes.at(i)).getChromosomename()<< " because it is shorter than "<< windowSize_ << endl;
		chromosomes_.erase(chromosomes_.begin()+smallChromosomes.at(i));

	}*/
}

int Genome::getNumberOfChromosomes() {
	return chromosomes_.size();
}

void Genome::addFromFasta(std::string const& genomeFastaFilename) {
	cout << "Reading the genome from "<< genomeFastaFilename <<endl;
	ifstream file;
	file.open(genomeFastaFilename.c_str());
	if(!file.is_open()) {
		file.clear();
		cerr << genomeFastaFilename <<"cannot be opened\n";
		exit (-1);
	}

	int chromosomeIndex = 0;
	int startOnChr = 0;
	char letter; // here we read the chromosome letter by letter. Can do it better!
	int count = 0;
	int countN = 0;
	string line;
    string text = "";
	string chromosome = "";

	while (!file.eof()) {
		file>>letter;
		if (letter == '>') {
			if (chromosome!="") { //save the previous segment if you want
				float notNpercentage = float(count-countN)/count;
				int left = startOnChr;
				int lastFragmentLength=text.length();
				startOnChr+=lastFragmentLength;
				chromosomes_.back().addElement(chromosome,notNpercentage,left,startOnChr-1,lastFragmentLength,text);
			}
			getline (file,chromosome);
			cout << "reading chromosome "<< chromosome << "\n";
			Chromosome newChromosome(windowSize_,chromosome);
			chromosomes_.push_back(newChromosome);
			//reset:
			countN = 0;
			count = 0;
			text="";
			startOnChr = 0;
		} else {
			countN += isN(letter);
			count ++;
			text+=letter;  //text.push_back(letter);
			if (count == windowSize_) {
				float notNpercentage = float(count-countN)/count;
				int left = startOnChr;
				startOnChr+=windowSize_;
				chromosomes_.back().addElement(chromosome,notNpercentage,left,startOnChr-1,windowSize_,text);
				//reset:
				countN = 0;
				count = 0;
				text="";
			}
		}
	}
	file.close();
	cout << "End-of-file reached.." << endl;
}

int Genome::getWindowSize(void) {
	return windowSize_;
}


//void Genome::fillCGprofile(std::string const& chrFolder) {
//	vector<GenomeElement>::iterator it;
//	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
//		it->fillCGprofile(chrFolder);
//	}
//}

Genome::~Genome(void)
{
	chromosomes_.clear();
}

int Genome::findChromForSV(int SVlength) {
	int n=chromosomes_.size();
	int selectedChr= rand() % n;
	cout << "Selecting randomly one out of "<< n<< " chromosomes"<<endl;
	if (chromosomes_[selectedChr].getLengthBins()*windowSize_ < SVlength*1.5 || chromosomes_[selectedChr].getLengthBins()<MINCHRLENFORSV) { //try to find a longer chromosome:
		selectedChr = -1;
		for (int trial = 0; trial < 15000; trial ++) {
			selectedChr= rand() % n;
			if (chromosomes_[selectedChr].getLengthBins()*windowSize_ < SVlength*1.5 || chromosomes_[selectedChr].getLengthBins()<MINCHRLENFORSV) {
				//cout << "Selected chromosome with "<< chromosomes_[selectedChr].getLengthBins()<< " bins" <<endl;
				selectedChr = -1;
			} else {
				break;
			}
		}
	}
	if (selectedChr==-1)
		cout <<"Could not find a chromosome in "<<15000<<" trials..";
	return selectedChr;
}

void Genome::addInversions(int number, int minLength, int maxLength){
	cout << "..Adding "<< number <<" inversions from "<<minLength<<" to " <<maxLength<< " bp\n";
	vector<int> insertSize (number);
	fillSizes(insertSize,number,minLength,maxLength);

	for (int i=0; i<number; i++) {
		//choose randomly a chromosome where to insert it
		int selectedChr = findChromForSV(insertSize[i]);
		if (selectedChr == -1) {
			cout << "..could not find a chromosome long enough to insert inversion of "<< insertSize[i] <<"bp" << endl;
			continue;
		}

		if (insertSize[i]<windowSize_) {
			chromosomes_[selectedChr].makeSmallInversion(insertSize[i]);
		} else {
			//invert whole elements:
			int numberOfElementsToDuplicate = insertSize[i] / windowSize_;
			chromosomes_[selectedChr].makeLargeInversion (numberOfElementsToDuplicate);
		}
	}
}

void Genome::addSmallInsertions(int number, int minLength, int maxLength) {
	cout << "..Adding "<< number <<" small insertions from "<<minLength<<" to " <<maxLength<< " bp\n";
	vector<int> insertSize (number);
	fillSizes(insertSize,number,minLength,maxLength);

	for (int i=0; i<number; i++) {
		//choose randomly a chromosome where to insert it
		int selectedChr = findChromForSV(insertSize[i]);
		if (selectedChr == -1) {
			cout << "..could not find a chromosome to insert "<< insertSize[i] <<"bp" << endl;
			continue;
		}
		chromosomes_[selectedChr].makeSmallInsertion(insertSize[i]);
	}
}
void Genome::addDuplicationsWithInversion(int number, int minLength, int maxLength) {
	cout << "..Adding "<< number <<" duplications with inversion from "<<minLength<<" to " <<maxLength<< " bp\n";
	vector<int> insertSize (number);
	fillSizes(insertSize,number,minLength,maxLength);
	for (int i=0; i<number; i++) {
		//choose randomly a chromosome where to insert it
		int selectedChr = findChromForSV(insertSize[i]);
		if (selectedChr == -1) {
			cout << "..could not find a chromosome to insert duplication with inversion of "<< insertSize[i] <<"bp" << endl;
			continue;
		}

		if (insertSize[i]<windowSize_) {
			chromosomes_[selectedChr].makeSmallDuplicationWithInversion(insertSize[i]);
		} else {
			//duplicate whole elements:
			int numberOfElementsToDuplicate = insertSize[i] / windowSize_;
			chromosomes_[selectedChr].makeLargeDuplicationWithInversion(numberOfElementsToDuplicate);
		}
	}

}
void Genome::addDuplications(int number, int minLength, int maxLength){
	cout << "..Adding "<< number <<" duplications from "<<minLength<<" to " <<maxLength<< " bp\n";
	vector<int> insertSize (number);
	fillSizes(insertSize,number,minLength,maxLength);

	for (int i=0; i<number; i++) {
		//choose randomly a chromosome where to insert it
		int selectedChr = findChromForSV(insertSize[i]);
		if (selectedChr == -1) {
			cout << "..could not find a chromosome to insert duplication of "<< insertSize[i] <<"bp" << endl;
			continue;
		}

		if (insertSize[i]<windowSize_) {
			chromosomes_[selectedChr].makeSmallDuplication(insertSize[i]);
		} else {
			//duplicate whole elements:
			int numberOfElementsToDuplicate = insertSize[i] / windowSize_;
			chromosomes_[selectedChr].makeLargeDuplication(numberOfElementsToDuplicate);
		}
	}
}

void Genome::addDeletions(int number, int minLength, int maxLength) {
	cout << "..Adding "<< number <<" deletions from "<<minLength<<" to " <<maxLength<< " bp\n";
	vector<int> insertSize (number);
	fillSizes(insertSize,number,minLength,maxLength);

	for (int i=0; i<number; i++) {
		//choose randomly a chromosome where to insert it
		int selectedChr = findChromForSV(insertSize[i]);
		if (selectedChr == -1) {
			cout << "..could not find a chromosome long enough to insert deletion of "<< insertSize[i] <<"bp" << endl;
			continue;
		}

		if (insertSize[i]<windowSize_) {
			chromosomes_[selectedChr].makeSmallSizeDeletion(insertSize[i]);
		} else {
			//delete whole elements:
			int numberOfElementsToDelete = insertSize[i] / windowSize_;
			chromosomes_[selectedChr].makeLargeDeletion(numberOfElementsToDelete);
		}
	}
}

void Genome::addLinkingInsertions(int number, int minLength, int maxLength, bool isInverted) {
	cout << "..Adding "<< number <<" linking insertions from "<<minLength<<" to " <<maxLength<< " bp\n";
	int n=chromosomes_.size();
	if (n==0) {
		cout << "..cannot create a linking insertion between two different chromosomes, since there is no chromosomes left" << endl;
		return;
	}

	vector<int> insertSize (number);
	fillSizes(insertSize,number,minLength,maxLength);

	for (int i=0; i<number; i++) {
		//choose randomly a chromosome where to insert it

		int selectedChr1= findChromForSV(1);
		int selectedChr2= findChromForSV(insertSize[i]);
		if (selectedChr2 == -1 || selectedChr1 == -1) {
			cout << "..could not find a chromosome to insert linking insertion of "<< insertSize[i] <<"bp" << endl;
			continue;
		}

		if (insertSize[i]<windowSize_) {
			GenomeElement * NewElement = chromosomes_[selectedChr2].getSmallFragment(insertSize[i]);
			if (NewElement!= NULL) {
				chromosomes_[selectedChr1].makeFragmentInsertion(NewElement,isInverted);
			}
		} else {
			//invert whole elements:
			int numberOfElementsToDuplicate = insertSize[i] / windowSize_;
			GenomeElement * NewElement = chromosomes_[selectedChr2].getLargeFragment(numberOfElementsToDuplicate);
			if (NewElement!= NULL) {
				chromosomes_[selectedChr1].makeFragmentInsertion(NewElement,isInverted);
			}
		}
	}
}

void Genome::addReInsertions(int number, int minLength, int maxLength, bool isInverted){
	cout << "..Adding "<< number <<" re-insertions from "<<minLength<<" to " <<maxLength<< " bp\n";
	int n=chromosomes_.size();
	if (n==0) {
		cout << "..cannot create a linking insertion between two different chromosomes, since there is no chromosomes left" << endl;
		return;
	}

	vector<int> insertSize (number);
	fillSizes(insertSize,number,minLength,maxLength);

	for (int i=0; i<number; i++) {
		//choose randomly a chromosome where to insert it

		int selectedChr1= findChromForSV(1);
		int selectedChr2= findChromForSV(insertSize[i]);
		if (selectedChr2 == -1 || selectedChr1==-1) {
			cout << "..could not find a chromosome to insert linking insertion of "<< insertSize[i] <<"bp" << endl;
			continue;
		}

		if (insertSize[i]<windowSize_) {
			GenomeElement * NewElement = chromosomes_[selectedChr2].extractSmallFragment(insertSize[i]);
			if (NewElement!= NULL) {
				chromosomes_[selectedChr1].makeFragmentInsertion(NewElement,isInverted);
			}
		} else {
			//insert whole elements:
			int numberOfElementsToDuplicate = insertSize[i] / windowSize_;
			GenomeElement * NewElement = chromosomes_[selectedChr2].extractLargeFragment(numberOfElementsToDuplicate);
			if (NewElement!= NULL) {
				chromosomes_[selectedChr1].makeFragmentInsertion(NewElement,isInverted);
			}
		}
	}

}

void Genome::addTranslocations(int number, bool isBalanced) {
	cout << "..Adding "<< number <<" translocations \n";
	int n=chromosomes_.size();
	if (n==1) {
		cout << "..cannot create a translocation between two different chromosomes, since there is only one chromosome left" << endl;
		return;
	}
	int selectedChr1= findChromForSV(1);
	int selectedChr2= findChromForSV(1);
	if (selectedChr1== -1 || selectedChr2==-1) {
		cout << "Cannot file chromosomes to maka a translocation"<<endl;
		return;
	}
	if (selectedChr1==selectedChr2) {
		bool whichOne = rand() % 2;
		if (selectedChr1!=0) {
			if (whichOne)
				selectedChr1--;
			else
				selectedChr2--;
		}else {
			if (whichOne)
				selectedChr1++;
			else
				selectedChr2++;
		}
	}

	bool isLeft1 = rand() % 2;
	bool isRight2 = rand() % 2;

	if (isBalanced) {
		Chromosome *newChrom1 =new Chromosome(windowSize_,chromosomes_[selectedChr1].getChromosomename()+"_"+chromosomes_[selectedChr2].getChromosomename());
		newChrom1->makeTransclocation(chromosomes_[selectedChr1],chromosomes_[selectedChr2],isLeft1,isRight2);
		Chromosome *newChrom2 =new Chromosome(windowSize_,chromosomes_[selectedChr1].getChromosomename()+"_"+chromosomes_[selectedChr2].getChromosomename());
		newChrom2->makeTransclocation(chromosomes_[selectedChr1],chromosomes_[selectedChr2],!isLeft1,!isRight2);
		chromosomes_.erase(chromosomes_.begin()+max(selectedChr1,selectedChr2));
		chromosomes_.erase(chromosomes_.begin()+min(selectedChr1,selectedChr2));
		chromosomes_.push_back(*newChrom1);
		chromosomes_.push_back(*newChrom2);
	}
	else {
		Chromosome *newChrom =new Chromosome(windowSize_,chromosomes_[selectedChr1].getChromosomename()+"_"+chromosomes_[selectedChr2].getChromosomename());
		newChrom->makeTransclocation(chromosomes_[selectedChr1],chromosomes_[selectedChr2],isLeft1,isRight2);
		bool isLossForFirst = rand() % 2;
		bool isLossForSecond = rand() % 2;
		if (isLossForFirst & isLossForSecond) {
			chromosomes_.erase(chromosomes_.begin()+max(selectedChr1,selectedChr2));
			chromosomes_.erase(chromosomes_.begin()+min(selectedChr1,selectedChr2));
		} else if (isLossForFirst) {
			chromosomes_.erase(chromosomes_.begin()+selectedChr1);
		} else if (isLossForSecond) {
			chromosomes_.erase(chromosomes_.begin()+selectedChr2);
		}
		chromosomes_.push_back(*newChrom);
	}

}

void Genome::deleteChromosomes(int number) {
	cout << "..Deleting "<< number<<" chromosomes\n";
	for (int i=0; i<number; i++) {
		//choose randomly a chromosome to delete
		if (chromosomes_.size()==0) {
			cout << "..your genome has lost all its chromosomes..."<< endl;
			return;
		}
		int selectedChr = findChromForSV(1);;
		chromosomes_.erase(chromosomes_.begin()+selectedChr);
	}
}

void Genome::duplicateChromosomes(int number){
	cout << "..Duplicating "<< number<<" chromosomes\n";
	if (chromosomes_.size()==0) {
			cout << "..your genome has lost all its chromosomes..."<< endl;
			return;
		}
	for (int i=0; i<number; i++) {
		//choose randomly a chromosome to duplicate
		int selectedChr = findChromForSV(1);;
		Chromosome *newChrom =new Chromosome(windowSize_,chromosomes_[selectedChr].getChromosomename());
		*newChrom =chromosomes_[selectedChr];
		chromosomes_.push_back(*newChrom);
	}
}

void Genome::duplicateGenome() {
	cout << "..Duplicating the genome\n";
	int n1=chromosomes_.size();
	if (n1==0) {
			cout << "..your genome has lost all its chromosomes..."<< endl;
			return;
	}
	for (int i=0; i<n1; i++) {
		int selectedChr = i;
		Chromosome *newChrom =new Chromosome(windowSize_,chromosomes_[selectedChr].getChromosomename());
		*newChrom =chromosomes_[selectedChr];
		chromosomes_.push_back(*newChrom);
	}
}

void Genome::printDonorGenome(std::string const& fastaFileName, std::string const& commentsFileName) {
	cout << "Printing the new tumor genome into "<< fastaFileName <<endl;
	std::ofstream file;
	file.open(fastaFileName.c_str());
	std::ofstream commentsFile;
	commentsFile.open(commentsFileName.c_str());
	int ID=0;
	for (std::vector<Chromosome>::iterator it=chromosomes_.begin() ; it != chromosomes_.end(); it++ ) {
		ID++;
		(*it).printDonorFasta(file,commentsFile,ID);
	}
	file.close();
	commentsFile.close();
}


void Genome::addCoAmplification(int number, int minLength, int maxLength, bool sameChr,int amplificationRate){
	cout << "..Adding "<< number <<" coamplifications from "<<minLength<<" to " <<maxLength<< " bp\n";
	int n=chromosomes_.size();
	if (n==0) {
		cout << "..cannot create an amplicon, since there is no chromosomes left" << endl;
		return;
	}
	string chrName = "Chromosome minute ";
	Chromosome *newChrom =new Chromosome(windowSize_,chrName);

	vector<int> insertSize (number);
	fillSizes(insertSize,number,minLength,maxLength);

	int selectedChr=1;
	if(sameChr) {
		selectedChr= findChromForSV(1);
		if (selectedChr==-1) {
			cout << "Cannot find chromosomes to create an amplicon"<<endl;
			return;
		}
		chrName+=chromosomes_[selectedChr].getChromosomename();
	}

	for (int i=0; i<number; i++) {
		if (!sameChr) {
			//choose randomly a chromosome and extract
			selectedChr= findChromForSV(1);
			if (selectedChr==-1) {
				cout << "Cannot find chromosomes to create an amplicon"<<endl;
				return;
			}
			chrName+=chromosomes_[selectedChr].getChromosomename();
		}
		if (selectedChr == -1) {
			cout << "..could not find a chromosome to extract an amplicon of "<< insertSize[i] <<"bp" << endl;
			return;
		}
		GenomeElement * NewElement;
		if (insertSize[i]<windowSize_) {
			NewElement = chromosomes_[selectedChr].getSmallFragment(insertSize[i]);
		} else {
			// several elements:
			int numberOfElementsToDuplicate = insertSize[i] / windowSize_;
			NewElement = chromosomes_[selectedChr].getLargeFragment(numberOfElementsToDuplicate);
		}

		if (NewElement!= NULL) {
				bool isInverted = rand() % 2;
				if(isInverted) {
					NewElement->invert();
				}
				newChrom->addElement(*NewElement);
		}
	}
	newChrom->amplifyMinute(amplificationRate);
	newChrom->setChromosomename(chrName);
	chromosomes_.push_back(*newChrom);
}
