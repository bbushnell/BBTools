package sketch;

import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Tools;
import structures.FloatList;
import tax.TaxNode;
import tax.TaxTree;

class ResultLineParser {

	ResultLineParser(int mode_, TaxTree tree_, ByteStreamWriter bswBad_, ArrayList<RecordSet> recordSets_, boolean keepText_){
		mode=mode_;
		tree=tree_;
		bswBad=bswBad_;
		recordSets=recordSets_;
		keepText=keepText_ || bswBad!=null;
		for(int i=0; i<AnalyzeSketchResults.taxLevels; i++){
			aniLists[i]=new FloatList();
			ssuLists[i]=new FloatList();
		}
	}

	void parse(byte[] line){
		if(keepText){text=line;}
		if(line[0]!='#'){
			if(mode==AnalyzeSketchResults.BBSKETCH_MODE){
				parseData(line);
			}else if(mode==AnalyzeSketchResults.MASH_MODE){
				parseDataMash(line);
			}else{
				assert(false) : "Bad mode: "+mode;
			}
		}else{
			parseHeader(line);
			if(bswBad!=null){bswBad.println(line);}
		}
	}

	private synchronized void parseHeader(byte[] line){
		ArrayList<byte[]> split=Tools.split(line, 0, (byte)'\t');
		for(int col=0; col<split.size(); col++){
			byte[] array=split.get(col);
			if(Tools.equals(array, "ANI") || Tools.equals(array, "AAI")){
				aniColumn=col;
			}else if(Tools.equals(array, "QTaxID")){
				qTaxIDColumn=col;
			}else if(Tools.equals(array, "RTaxID")){
				rTaxIDColumn=col;
			}else if(Tools.equals(array, "SSU")){
				ssuColumn=col;
			}else if(Tools.equals(array, "CALevel")){
				caLevelColumn=col;
			}

			else if(Tools.equals(array, "QSize")){
				qSizeColumn=col;
			}else if(Tools.equals(array, "RefSize") || Tools.equals(array, "RSize")){
				rSizeColumn=col;
			}else if(Tools.equals(array, "QBases")){
				qBasesColumn=col;
			}else if(Tools.equals(array, "RBases")){
				rBasesColumn=col;
			}
		}
	}

	private void parseData(byte[] line){
		ArrayList<byte[]> split=Tools.split(line, 0, (byte)'\t');
		qTaxID=Parse.parseInt(split.get(qTaxIDColumn), 0);
		rTaxID=Parse.parseInt(split.get(rTaxIDColumn), 0);
		qBases=Parse.parseLong(split.get(qBasesColumn), 0);
		rBases=Parse.parseLong(split.get(rBasesColumn), 0);
		qSize=Parse.parseLong(split.get(qSizeColumn), 0);
		rSize=Parse.parseLong(split.get(rSizeColumn), 0);
		ani=Parse.parseDouble(split.get(aniColumn), 0);
		byte[] ssuArray=split.get(ssuColumn);
		ssu=ssuArray.length==1 && ssuArray[0]=='.' ? -1 : Parse.parseDouble(ssuArray, 0);
		taxLevelExtended=TaxTree.stringToLevelExtended(new String(split.get(caLevelColumn)));
		if(taxLevelExtended<0) {
			System.err.println(new String(split.get(caLevelColumn)));
			taxLevelExtended=0;
		}
		processed=false;
	}
	
	private TaxNode getTaxNode(String fname){
		String name=ReadWrite.stripToCore(fname);
		if(name.startsWith("tid_")){
			int idx2=fname.indexOf('_', 4);
			int x=Parse.parseInt(fname, 4, idx2);
			return x>0 ? tree.getNode(x) : null;
			//name=name.substring(idx2+1); //This would allow fall-through to name parsing
		}
		try {
			return tree.getNodeByName(name);
		} catch (Throwable e) {
			return null;
		}
	}

	private void parseDataMash(byte[] line){
		///dev/shm/tid_123_Zymomonas_mobilis.fna.gz	/dev/shm/tid_456_bacterium_endosymbiont_of_Bathymodiolus_sp._5_South.fna.gz	0.43859	0.00515848	1/20000

		String[] split=new String(line).split("\t");

		String fraction=split[split.length-1];
		int numerator=Integer.parseInt(fraction.split("/")[0]);
		if(numerator<MIN_HITS){return;}
		int denominator=Integer.parseInt(fraction.split("/")[1]);

		//The default ordering is reversed since mash output is ordered first by ref, then query
		//The normal ordering (as below) requires a linux sort
		{
			TaxNode qNode=getTaxNode(split[0]);
			TaxNode rNode=getTaxNode(split[1]);

			if(qNode==null || rNode==null){return;}
			qTaxID=qNode.id;
			rTaxID=rNode.id;
			TaxNode ancestor=tree.commonAncestor(qNode, rNode);
			taxLevelExtended=ancestor.levelExtended;
		}
		
		ani=numerator/(float)denominator;
		ssu=-1;
		if(taxLevelExtended<0){taxLevelExtended=0;}
		processed=false;
	}

	//Returns a complete set when a new set is started
	RecordSet processData(HashMap<Long, Float> map, boolean saveRecord){
		RecordSet old=null;
		if(processed){return null;}
		levelAniSums[taxLevelExtended]+=ani;
		levelCounts[taxLevelExtended]++;
		aniLists[taxLevelExtended].add((float)ani);

		if(ssu>0){
			levelSSUSums[taxLevelExtended]+=ssu;
			levelCountsSSU[taxLevelExtended]++;
			ssuLists[taxLevelExtended].add((float)ssu);
		}
		if(map!=null){
			long key=(((long)qTaxID)<<32)|rTaxID;
			map.put(key, (float)ani);
		}
		if(saveRecord){
			if(currentSet==null || currentSet.qID!=qTaxID){
				old=currentSet;
				currentSet=new RecordSet(qTaxID);
				if(recordSets!=null){
					recordSets.add(currentSet);
				}
			}
			currentSet.records.add(new Record(this));
		}
		processed=true;
		return old;
	}

	/*--------------------------------------------------------------*/

	//		final static int taxLevels=TaxTree.numTaxaNamesExtended;
	final long[] levelCounts=new long[AnalyzeSketchResults.taxLevels];
	final long[] levelCountsSSU=new long[AnalyzeSketchResults.taxLevels];

	final double[] levelAniSums=new double[AnalyzeSketchResults.taxLevels];
	final double[] levelSSUSums=new double[AnalyzeSketchResults.taxLevels];

	final FloatList[] aniLists=new FloatList[AnalyzeSketchResults.taxLevels];
	final FloatList[] ssuLists=new FloatList[AnalyzeSketchResults.taxLevels];

	final ArrayList<RecordSet> recordSets;

	final int mode;
	final TaxTree tree;
	final ByteStreamWriter bswBad;

	int qTaxID=-1;
	int rTaxID=-1;
	long qBases;
	long rBases;
	long qSize;
	long rSize;
	double ani=-1;
	double ssu=-1;
	int taxLevelExtended=-1;
	boolean processed=true;
	RecordSet currentSet=null;
	final boolean keepText;

	byte[] text=null;

	private static int qTaxIDColumn=7;
	private static int rTaxIDColumn=8;
	private static int qSizeColumn=3;
	private static int rSizeColumn=4;
	private static int qBasesColumn=5;
	private static int rBasesColumn=6;
	private static int aniColumn=2;
	private static int ssuColumn=11;
	private static int caLevelColumn=12;

	static int MIN_HITS=3;
	
}