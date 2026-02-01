package sketch;

import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Tools;
import structures.FloatList;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Parses tabular result lines from genomic comparison sketches, extracting taxonomic
 * and similarity metrics. Supports both BBSketch and MASH result formats with
 * comprehensive data extraction capabilities including ANI, SSU, and taxonomic information.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
class ResultLineParser {

	/**
	 * Constructs a parser for genomic sketch comparison results.
	 * Initializes parsing mode, taxonomy tree, and data collection structures.
	 *
	 * @param mode_ Parsing mode (BBSKETCH_MODE or MASH_MODE)
	 * @param tree_ Taxonomic tree for resolving taxonomic relationships
	 * @param bswBad_ Stream writer for outputting bad/header lines (may be null)
	 * @param recordSets_ Collection to store parsed record sets (may be null)
	 * @param keepText_ Whether to preserve original text lines
	 */
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

	/**
	 * Parses a single result line based on the configured mode.
	 * Handles both data lines and header lines, routing to appropriate parsing methods.
	 * @param line Raw byte array containing the line to parse
	 */
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

	/**
	 * Parses header line to identify column positions for data extraction.
	 * Maps column names to indices for ANI, QTaxID, RTaxID, SSU, CALevel, and size columns.
	 * Thread-safe to handle concurrent access during multi-threaded processing.
	 * @param line Header line containing tab-separated column names
	 */
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

	/**
	 * Parses BBSketch format data line extracting all comparison metrics.
	 * Extracts taxonomic IDs, sequence sizes, ANI values, SSU scores, and taxonomic levels.
	 * @param line Tab-separated data line in BBSketch format
	 */
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
	
	/**
	 * Extracts taxonomic node from filename, supporting both TaxID and name-based lookup.
	 * Handles filename formats with "tid_" prefix for direct taxonomic ID extraction.
	 * @param fname Filename containing taxonomic information
	 * @return TaxNode corresponding to the filename, or null if not found
	 */
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

	/**
	 * Parses MASH format data line extracting comparison metrics.
	 * Processes filename-based taxonomic identification and fraction-based similarity scoring.
	 * Filters results based on minimum hit threshold.
	 * @param line Tab-separated data line in MASH format
	 */
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
	/**
	 * Processes parsed data by updating statistics and optionally saving records.
	 * Accumulates ANI and SSU statistics by taxonomic level and manages record sets.
	 *
	 * @param map Optional hash map to store query-reference ANI pairs
	 * @param saveRecord Whether to create and store Record objects
	 * @return Previous RecordSet if a new query ID was encountered, null otherwise
	 */
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
	/** Count of comparisons at each taxonomic level */
	final long[] levelCounts=new long[AnalyzeSketchResults.taxLevels];
	/** Count of SSU comparisons at each taxonomic level */
	final long[] levelCountsSSU=new long[AnalyzeSketchResults.taxLevels];

	/** Sum of ANI values at each taxonomic level for averaging */
	final double[] levelAniSums=new double[AnalyzeSketchResults.taxLevels];
	/** Sum of SSU values at each taxonomic level for averaging */
	final double[] levelSSUSums=new double[AnalyzeSketchResults.taxLevels];

	/** Lists of individual ANI values at each taxonomic level */
	final FloatList[] aniLists=new FloatList[AnalyzeSketchResults.taxLevels];
	/** Lists of individual SSU values at each taxonomic level */
	final FloatList[] ssuLists=new FloatList[AnalyzeSketchResults.taxLevels];

	/** Collection of record sets for storing parsed comparison records */
	final ArrayList<RecordSet> recordSets;

	/** Parsing mode indicating format type (BBSKETCH_MODE or MASH_MODE) */
	final int mode;
	/** Taxonomic tree for resolving taxonomic relationships and ancestry */
	final TaxTree tree;
	/** Stream writer for outputting header lines or bad data (may be null) */
	final ByteStreamWriter bswBad;

	/** Query taxonomic ID from current parsed line */
	int qTaxID=-1;
	/** Reference taxonomic ID from current parsed line */
	int rTaxID=-1;
	/** Number of bases in query sequence */
	long qBases;
	/** Number of bases in reference sequence */
	long rBases;
	/** Size of query sketch in k-mers */
	long qSize;
	/** Size of reference sketch in k-mers */
	long rSize;
	/** Average Nucleotide Identity value from current comparison */
	double ani=-1;
	/** Small Subunit (16S/18S) rRNA similarity score (-1 if not available) */
	double ssu=-1;
	/** Extended taxonomic level of common ancestor between query and reference */
	int taxLevelExtended=-1;
	/** Flag indicating whether current data has been processed */
	boolean processed=true;
	/** Current record set being populated with comparison records */
	RecordSet currentSet=null;
	/** Whether to preserve original text of parsed lines */
	final boolean keepText;

	/** Original text of current line if keepText is enabled */
	byte[] text=null;

	/** Column index for query taxonomic ID in result files */
	private static int qTaxIDColumn=7;
	/** Column index for reference taxonomic ID in result files */
	private static int rTaxIDColumn=8;
	/** Column index for query sketch size in result files */
	private static int qSizeColumn=3;
	/** Column index for reference sketch size in result files */
	private static int rSizeColumn=4;
	/** Column index for query sequence base count in result files */
	private static int qBasesColumn=5;
	/** Column index for reference sequence base count in result files */
	private static int rBasesColumn=6;
	/** Column index for ANI (Average Nucleotide Identity) values in result files */
	private static int aniColumn=2;
	/**
	 * Column index for SSU (Small Subunit rRNA) similarity values in result files
	 */
	private static int ssuColumn=11;
	/** Column index for common ancestor taxonomic level in result files */
	private static int caLevelColumn=12;

	/** Minimum number of k-mer hits required for MASH format processing */
	static int MIN_HITS=3;
	
}