#include "Utils.h"
#include "BiopolymerClass.h"

typedef char TChar;                             // character type
typedef seqan::String<TChar> TSequence;                // sequence type
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;



struct ThreadingPartner {
    BiopolymerClass biopolymerClass;
    vector <ResidueID> includedResidues;
    ResidueID startResidue;
    ResidueID endResidue;
    TSequence sequence;
};


class ThreadingStruct {
   private:
       ThreadingPartner threadingPartners[2];
       seqan::AlignmentStats alignmentStats;
       TAlign align;
       bool alignHasBeenComputed;
       
    public:
        //double gapPenalty;
        // This method is not needed. just use updThreadingPartner.
        //void setThreadingPartner (ThreadingPartner myThreadingPartner, int index){threadingPartners[index] =  myThreadingPartner;}
        // In homologyScanner, per convention partner 0 is the homologJob, partner 1 is the PrimaryJob
        ThreadingPartner & updThreadingPartner (int index){return threadingPartners[index];}
        ThreadingPartner  getThreadingPartner (int index) const {return threadingPartners[index];}
        std::string getChain(int index){return threadingPartners[index].biopolymerClass.getChainID();};
        //void print(){}
        void setDefaultStartEndResidues(){
            updThreadingPartner(0).startResidue = updThreadingPartner(0).biopolymerClass.getFirstResidueID();
            updThreadingPartner(0).  endResidue = updThreadingPartner(0).biopolymerClass. getLastResidueID();
            updThreadingPartner(1).startResidue = updThreadingPartner(1).biopolymerClass.getFirstResidueID();
            updThreadingPartner(1).  endResidue = updThreadingPartner(1).biopolymerClass. getLastResidueID();
        }
        double forceConstant;
        bool backboneOnly;
        bool isGapped; //if False, then alignment is being provided explicitly. if True, precise alignment will be determined by MMB/SeqAn
        bool deadLengthIsFractionOfInitialLength; //If True, then dead length of each spring will be set to deadLengthFraction * <initial spring extension>. It makes sense that 1 > deadLengthFraction > 0.
        double deadLength; // This is an absolute dead length for the alignment springs. For default homology modeling behavior, should be 0.
        double deadLengthFraction;
        seqan::AlignmentStats getAlignmentStats(){return alignmentStats;}
        void   setAlignmentStats(seqan::AlignmentStats myAlignmentStats){ alignmentStats = myAlignmentStats;}

        // align, alignmentStats should be empty, undefined, or garbage unless and until this method is called:
        void setLongSequences(){
            threadingPartners[0].sequence = (threadingPartners[0].biopolymerClass.getSubSequence(threadingPartners[0].startResidue , threadingPartners[0].endResidue )).c_str();
            threadingPartners[1].sequence = (threadingPartners[1].biopolymerClass.getSubSequence(threadingPartners[1].startResidue , threadingPartners[1].endResidue )).c_str();
            computeAlign(); // Just to make sure align is in sync with the new sequences
        }

        void sortIncludedResidues(){
            //threadingPartners[0].includedResidues = 
            threadingPartners[0].biopolymerClass.sort(threadingPartners[0].includedResidues);
            //threadingPartners[1].includedResidues = 
            threadingPartners[1].biopolymerClass.sort(threadingPartners[1].includedResidues);
        }

	bool hasResidue( ResidueID residue , int biopolymerIndex ){
	    for (int i = 0; i < threadingPartners[biopolymerIndex].includedResidues.size() ;  i++){
		if  (  threadingPartners[biopolymerIndex].includedResidues[i] == residue) return 1;
	    }
	    return 0;
	}

	void supplementIncludedResidues(int fromBiopolymer, int toBiopolymer){
	    for (int i = 0; i < threadingPartners[fromBiopolymer].includedResidues.size() ;  i++){
		ResidueID correspondingResidue;
		if (!( getCorrespondingResidue(threadingPartners[fromBiopolymer].includedResidues[i], correspondingResidue, fromBiopolymer, toBiopolymer) )) // returns 0 if successful
		{
		    if (!(hasResidue(  correspondingResidue , toBiopolymer ) )){
			threadingPartners[toBiopolymer].includedResidues.push_back(correspondingResidue);
		    }  // of if
		} // of if
	    }  // of for
	} // of method

	void supplementIncludedResidues(){
	    supplementIncludedResidues(0,1);
	    supplementIncludedResidues(1,0);
	} // of method

	void setShortSequences(){
	    sortIncludedResidues();
	    threadingPartners[0].sequence = (threadingPartners[0].biopolymerClass.getSequence(threadingPartners[0].includedResidues)).c_str();
	    threadingPartners[1].sequence = (threadingPartners[1].biopolymerClass.getSequence(threadingPartners[1].includedResidues)).c_str();
	    computeAlign(); // Just to make sure align is in sync with the new sequences
	}

        TAlign computeAlign(){
	    seqan::resize(rows(align), 2);
	    assignSource(row(align,0),threadingPartners[0].sequence);
	    assignSource(row(align,1),threadingPartners[1].sequence);
	    seqan::Blosum62 scoringScheme(-1, -12);
	    int score = globalAlignment(align,scoringScheme ); // ..signature:Score<TValue, Simple>(match, mismatch, gap [, gap_open])
	    //"62" means matrix was constructed with max 62% seq. identity alignment.
	    std::cout <<__FILE__<<":"<<__LINE__<< " SeqAn sequence alignment follows: "  << ::std::endl;
	    std::cout <<__FILE__<<":"<<__LINE__<< "Score: " << score << ::std::endl;
	    std::cout <<__FILE__<<":"<<__LINE__<< align << ::std::endl;
	    computeAlignmentStats(alignmentStats, align, scoringScheme);
	    printAlignmentStats();
	    alignHasBeenComputed = 1;
            return align;
	}

	// return 0 for success, 1 for failure
	bool getCorrespondingResidue(const ResidueID queryResidue, ResidueID & correspondingResidue, const int queryBiopolymerIndex, const int correspondingBiopolymerIndex){
	    if (!( queryBiopolymerIndex ^ correspondingBiopolymerIndex )){std::cout <<__FILE__<<":"<<__LINE__<< " queryBiopolymerIndex and  correspondingBiopolymerIndex must be set to 0,1 or 1,0. You have   "<< queryBiopolymerIndex <<" , " << correspondingBiopolymerIndex <<std::endl; exit(1);}
	    if (( queryBiopolymerIndex > 1) ||  ( queryBiopolymerIndex < 0) ){std::cout <<__FILE__<<":"<<__LINE__<< " queryBiopolymerIndex =  "<< queryBiopolymerIndex <<" is out of range. Expected 0 or 1. "  <<std::endl; exit(1);}
	    if (( correspondingBiopolymerIndex > 1) ||  ( correspondingBiopolymerIndex < 0) ){std::cout <<__FILE__<<":"<<__LINE__<< " correspondingBiopolymerIndex =  "<< correspondingBiopolymerIndex <<" is out of range. Expected 0 or 1. "  <<std::endl; exit(1);}
	    if (!(alignHasBeenComputed)) {std::cout <<__FILE__<<":"<<__LINE__<< " alignHasBeenComputed = "<< alignHasBeenComputed <<std::endl; exit(1);}
	    int querySourceIndex = threadingPartners[queryBiopolymerIndex].biopolymerClass.getResidueIndex(queryResidue) -  threadingPartners[queryBiopolymerIndex].biopolymerClass.getResidueIndex(threadingPartners[queryBiopolymerIndex].startResidue);
	     int viewIndex = toViewPosition(row(align, queryBiopolymerIndex),querySourceIndex);
	     int correspondingSourceIndex =  toSourcePosition(row(align,correspondingBiopolymerIndex ) ,viewIndex);
	     if (String(seqan::row (align,correspondingBiopolymerIndex )[viewIndex]) == ("-")   ){
		 //std::cout <<__FILE__<<":"<<__LINE__<< " queryResidue = "<<queryResidue.outString()<< " of queryBiopolymerIndex "<<queryBiopolymerIndex<< " with chain ID "<<   threadingPartners[queryBiopolymerIndex].biopolymerClass.getChainID() << " is an insertion with respect to correspondingBiopolymerIndex. Unable to set correspondingResidue "  <<std::endl;
		 return 1;
	     }
	     correspondingResidue =  threadingPartners[correspondingBiopolymerIndex].biopolymerClass.sum( threadingPartners[correspondingBiopolymerIndex].startResidue, correspondingSourceIndex);
	     return 0;
	 }; // of method

	 const void printAlignmentStats(){
	     std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
	     std::cout << align
		   << "score:               " << alignmentStats.alignmentScore << "\n"
		   << "gap opens:           " << alignmentStats.numGapOpens << "\n"
		   << "gap extensions:      " << alignmentStats.numGapExtensions << "\n"
		   << "num insertions:      " << alignmentStats.numInsertions << "\n"
		   << "num deletions:       " << alignmentStats.numDeletions << "\n"
		   << "num matches:         " << alignmentStats.numMatches << "\n"
		   << "num mismatches:      " << alignmentStats.numMismatches << "\n"
		   << "num positive scores: " << alignmentStats.numPositiveScores << "\n"
		   << "num negative scores: " << alignmentStats.numNegativeScores << "\n"
		   << "percent similarity:  " << alignmentStats.alignmentSimilarity << "\n"
		   << "alignment length  :  " << alignmentStats.alignmentLength     << "\n"
		   << "percent identity:    " << alignmentStats.alignmentIdentity << "\n\n\n";
	     std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
	 }


	ThreadingStruct(BiopolymerClass  biopolymerClass0, ResidueID resStart0, ResidueID resEnd0,
		     BiopolymerClass  biopolymerClass1, ResidueID resStart1, ResidueID resEnd1,
		     double forceConstant, bool backboneOnly
		     ) //:

	    //residueStart0(resStart0), residueEnd0(resEnd0),
			 //residueStart1(resStart1), residueEnd1(resEnd1),
			 //forceConstant(forceConstant), backboneOnly(backboneOnly)
	    {   
		updThreadingPartner(0).biopolymerClass = biopolymerClass0; 
		updThreadingPartner(1).biopolymerClass = biopolymerClass1; 
		updThreadingPartner(0).includedResidues.clear();
		//chain1ResiduesIncluded.clear(); 
		alignHasBeenComputed = 0; 
		threadingPartners[0].sequence = ""; threadingPartners[1].sequence = "";}

		// This constructor has to change because we are getting rido of chainID's. Also, usage in MMB has to change.
	ThreadingStruct() 
		     {
	        alignHasBeenComputed = 0;
	        threadingPartners[0].includedResidues.clear(); threadingPartners[1].includedResidues.clear();
                threadingPartners[0].sequence = ""; threadingPartners[1].sequence = "";}
}; // of class


