package bin;

import aligner.SingleStateAlignerFlat2;
import dna.Data;
import prok.CallGenes;
import prok.GeneCaller;
import prok.GeneModel;
import prok.GeneModelParser;

public class GeneTools {
	
	static synchronized GeneCaller makeGeneCaller() {
		if(pgm==null) {loadPGM();}
		return CallGenes.makeGeneCaller(pgm);
	}
	
	static synchronized void loadPGM() {
		if(pgm!=null) {return;}
		if(pgmFile==null){pgmFile=Data.findPath("?model.pgm");}
		
		if(!quiet) {System.err.println("Loading "+pgmFile);}
		CallGenes.call16S=CallGenes.call18S=true;
		CallGenes.loadLongKmers();
		CallGenes.loadConsensusSequenceFromFile(true, true);
		pgm=GeneModelParser.loadModel(pgmFile);
		gCaller=(pgm==null && gCaller!=null ? null : CallGenes.makeGeneCaller(pgm));
	}
	
	static synchronized void setMode(boolean r16, boolean r18, boolean r5, boolean r23, boolean trna, boolean cds) {
		GeneCaller.call16S=r16;
		GeneCaller.call18S=r18;
		GeneCaller.call5S=r5;
		GeneCaller.call23S=r23;
		GeneCaller.calltRNA=trna;
		GeneCaller.callCDS=cds;
	}
	
	static String pgmFile=Data.findPath("?model.pgm");
	static GeneModel pgm;
	static GeneCaller gCaller;
	static boolean quiet=false;
	
}
