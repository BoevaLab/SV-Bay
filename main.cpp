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

// main.cpp

#include "tgsim.h"
#include "Genome.h"

using namespace std ;
int verbose = false;

int main(int argc, char *argv[]) {
	//check arguments
    cout << "Tumor Genome Simulator (TGSim), version 1.0\n";
    cout << "TGSim adds large structural variants to a normal genome represented by a fasta file.\nCopyright (C) 2014 Institut Curie <http://www.curie.fr>\n";
	if (argc < 3) {
		cerr << "\n\tPlease specify a config file:\n\n./tgsim -conf config.txt\n";
		return 0;
	}
	if (argc == 3) {
		if (
			(strcmp(argv[1], "-conf") != 0)&&
			(strcmp(argv[1], "--conf") != 0)
			) {
            cerr << "\n\tPlease specify a config file\n\n";
			return 0;
		}
	}
	//check if config exists
    ifstream ifile(argv[2]);
    if (!ifile) {
      cerr << "\n\tCould not find your config file.. Please, check the existance of "<<argv[2] <<"\n\n";
      return -1;
    }

	//read config file
	ConfigFile cf(argv[2]);

    int window;

	try {
		window = int(cf.Value("general","window"));
        cout << "..windows size was set to "<< window<< "\n";

	} catch (const char * error) {
        cerr << "you need to specify window size\n";
        exit(0);
	}
	string outputFastaName = "donor.fa";

	string outputDir = ".";
	try {
		outputDir = std::string(cf.Value("general","outputFastaName"));
	} catch (const char * error) {
		//Do nothing, outputFastaName will be "donor.fa"
	}



	try {
		outputDir = std::string(cf.Value("general","outputDir"));
        if ( access( outputDir.c_str(), 0 ) == 0 )
        {
            struct stat status;
            stat( outputDir.c_str(), &status );

            if ( status.st_mode & S_IFDIR )
            {
                cout << "..Output directory:\t" << outputDir << "\n";
            }
            else
            {
                cerr << "Error: The path you entered for 'outputDir': "<< outputDir <<" is a file. It shoud be a directory" << endl;
                exit(-1);
            }
        }
        else
        {
            cerr << "Error: Path "<<outputDir<<" doesn't exist." << endl;
            exit(-1);
        }


	} catch (const char * error) {
		//Do nothing, outputDir will be "."
	}
	char sep = '/';
	#ifdef _WIN32
	   sep = '\\';
	#endif
	string commentFile = outputDir + sep + "comments_for_" + outputFastaName + ".txt";
	outputFastaName = outputDir + sep + outputFastaName;
    string fastaSeq = "";

    try {
            fastaSeq = std::string(cf.Value("general","genomeFasta"));

            if ( access( fastaSeq.c_str(), 0 ) == 0 )
            {
                struct stat status;
                stat( fastaSeq.c_str(), &status );

                if ( status.st_mode & S_IREAD )
                {
                    cout << "..File containing chromosome sequences:\t" << fastaSeq << "\n";
                }
                else
                {
                    cerr << "Error: The path you entered for 'fastaSeq': "<< fastaSeq <<" is a not file." << endl;
                    exit(-1);
                }
            }
            else
            {
                cerr << "Error: Path "<<fastaSeq<<" doesn't exist." << endl;
                exit(-1);
            }
    } catch (const char * error) {
            cerr << "..you need to provide a file with fasta sequences. Unable to proceed\n";
            exit(-1);
    }

	 string fastaSeq2 = "";

	 try {
            fastaSeq2 = std::string(cf.Value("general","genomeFasta2"));

            if ( access( fastaSeq2.c_str(), 0 ) == 0 )
            {
                struct stat status;
                stat( fastaSeq2.c_str(), &status );

                if ( status.st_mode & S_IREAD )
                {
                    cout << "..Second file containing chromosome sequences:\t" << fastaSeq2 << "\n";
                }
                else
                {
                    cerr << "Error: The path you entered for 'fastaSeq': "<< fastaSeq2 <<" is a not file." << endl;
                    exit(-1);
                }
            }
            else
            {
                cerr << "Error: Path "<<fastaSeq2<<" doesn't exist." << endl;
                exit(-1);
            }
    } catch (const char * error) {
           //do nothing
    }

    //READ SAMPLE DATA:

	Genome genomeInstance (window,fastaSeq);
	if (fastaSeq2!="")
		genomeInstance.addFromFasta(fastaSeq2);

	cout << "The genome consist of "<< genomeInstance.getNumberOfChromosomes()<< " chromosomes\n";

	//ADD SVS:
	genomeInstance.addInversions(3,6000,750000);
	genomeInstance.addDuplications(2,3000,10000);
	genomeInstance.addDeletions(2,10000,700000);
	genomeInstance.duplicateGenome();
	genomeInstance.deleteChromosomes(8);
	genomeInstance.duplicateChromosomes(2);
	genomeInstance.addTranslocations(1,BALANCED);
	genomeInstance.addTranslocations(4,UNBALANCED);
	genomeInstance.addDuplications(5,3000,1200000);
	genomeInstance.addDeletions(3,100000,3000000);
	genomeInstance.addDeletions(3,3000,100000); //small size after the large size
	genomeInstance.addLinkingInsertions(3,4000,40000,DIRECT);
	genomeInstance.addLinkingInsertions(2,10000,30000,INVERTED);
	genomeInstance.addReInsertions(2,3000,80000,DIRECT);
	genomeInstance.addReInsertions(2,4000,60000,INVERTED);
	genomeInstance.addDuplicationsWithInversion(2,3000,40000);
	genomeInstance.addInversions(4,3000,200000);
	genomeInstance.addSmallInsertions(5,1000,2000);
	genomeInstance.duplicateChromosomes(1);
	genomeInstance.deleteChromosomes(1);
	genomeInstance.addTranslocations(1,UNBALANCED);
	genomeInstance.duplicateChromosomes(1);

	//add amplicons:
	genomeInstance.addCoAmplification(4,60000,80000,SAMECHR,20);
	genomeInstance.addCoAmplification(3,60000,80000,DIFFCHR,50);

//whole chromosome duplications and deletions
	genomeInstance.printDonorGenome(outputFastaName,commentFile);
	return 0;
}
