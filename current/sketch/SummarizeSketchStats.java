package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Colors;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Tools;
import tax.TaxNode;
import tax.TaxTree;

/**
 * @author Brian Bushnell
 * @date June 28, 2017
 *
 */
public class SummarizeSketchStats {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Create a new SummarizeSketchStats instance
		SummarizeSketchStats x=new SummarizeSketchStats(args);
		
		///And run it
		x.summarize();
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs SummarizeSketchStats with argument parsing.
	 * Parses command line arguments for input files, taxonomic filters,
	 * and output configuration options.
	 * @param args Command line arguments
	 */
	public SummarizeSketchStats(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		ArrayList<String> names=new ArrayList<String>();
		String taxTreeFile=null;
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("printtotal") || a.equals("pt")){
				printTotal=Parse.parseBoolean(b);
			}else if(a.equals("ignoresametaxa")){
				ignoreSameTaxa=Parse.parseBoolean(b);
			}else if(a.equals("ignoresamebarcode") || a.equals("ignoresameindex")){
				ignoreSameBarcode=Parse.parseBoolean(b);
			}else if(a.equals("ignoresamelocation") || a.equals("ignoresameloc")){
				ignoreSameLocation=Parse.parseBoolean(b);
			}else if(a.equals("usetotal") || a.equals("totaldenominator") || a.equals("totald") || a.equals("td")){
				totalDenominator=Parse.parseBoolean(b);
			}
			
			else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("level") || a.equals("lv") || a.equals("taxlevel") || a.equals("tl") || a.equals("minlevel")){
				taxLevel=TaxTree.parseLevel(b);
				if(taxLevel>=0){
					taxLevel=TaxTree.levelToExtended(taxLevel);
				}
			}else if(a.equalsIgnoreCase("unique") || a.equalsIgnoreCase("uniquehits")){
				uniqueHitsForSecond=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("header") || a.equalsIgnoreCase("printheader")){
				printHeader=Parse.parseBoolean(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(!arg.contains("=")){
				String[] x=(new File(arg).exists() ? new String[] {arg} : arg.split(","));
				for(String x2 : x){names.add(x2);}
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		
		{//Process parser fields
			out=(parser.out1==null ? "stdout" : parser.out1);
			if(parser.in1!=null){
				String[] x=(new File(parser.in1).exists() ? new String[] {parser.in1} : parser.in1.split(","));
				for(String x2 : x){names.add(x2);}
			}
		}

		in=new ArrayList<String>();
		for(String s : names){
			Tools.getFileOrFiles(s, in, false, false, false, true);
		}
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile);}
	}
	
	/** Loads taxonomic tree from specified file.
	 * @param taxTreeFile Path to taxonomic tree file, or null to skip loading */
	void setTaxtree(String taxTreeFile){
		if(taxTreeFile==null){
			return;
		}
		tree=TaxTree.loadTaxTree(taxTreeFile, outstream, false, false);
	}
	
	/**
	 * Processes all input files and generates summary output.
	 * Reads sketch result files, creates summaries, and writes
	 * consolidated output with optional header.
	 */
	public void summarize(){
		ArrayList<SketchResultsSummary> list=new ArrayList<SketchResultsSummary>();
		for(String fname : in){
			ArrayList<SketchResultsSummary> ssl=summarize(fname);
			list.addAll(ssl);
		}
		
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		if(printHeader){tsw.print(header());}
//		if(printTotal){
//			tsw.println(total.toString());
//		}
		for(SketchResultsSummary ss : list){
			tsw.print(ss.toString());
		}
		tsw.poisonAndWait();
	}
	
//	Query: Troseus_1X_k55.fa	Seqs: 121 	Bases: 2410606	gSize: 2368581	SketchLen: 8923
//	WKID	KID	ANI	Complt	Contam	Matches	Unique	noHit	TaxID	gSize	gSeqs	taxName
//	99.89%	50.73%	100.00%	50.77%	0.02%	5683	5683	5	0	4719674	1	.	Troseus
	
	/**
	 * Parses a single sketch results file into summary objects.
	 * Processes file line by line, extracting query headers and
	 * associated result lines into SketchResultsSummary objects.
	 * @param fname Path to sketch results file
	 * @return List of summary objects from the file
	 */
	private ArrayList<SketchResultsSummary> summarize(String fname){
		TextFile tf=new TextFile(fname);
		ArrayList<SketchResultsSummary> list=new ArrayList<SketchResultsSummary>();
		SketchResultsSummary current=null;
		
		final String format="WKID	KID	ANI	Complt	Contam	Matches	Unique	noHit	TaxID	gSize	gSeqs	taxName";
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("Query:")){
				if(current!=null){list.add(current);}
				current=new SketchResultsSummary(line);
			}else if(line.startsWith("WKID")){
				assert(line.equals(format)) :
					"Format should be:\n"+format;
			}else if(line.length()>0){
				assert(current!=null) : "No Query Header for line "+line;
				current.add(line);
			}
		}
		if(current!=null){list.add(current);}
		tf.close();
		return list;
	}
	
	/** Generates tab-delimited header line for output format.
	 * @return Header string with column names for summary output */
	public static String header(){
		StringBuilder sb=new StringBuilder();
		
		sb.append("#query");

		sb.append('\t').append("seqs");
		sb.append('\t').append("bases");
		sb.append('\t').append("gSize");
		sb.append('\t').append("sketchLen");
		
		sb.append('\t').append("primaryHits");
		sb.append('\t').append("primaryUnique");
		sb.append('\t').append("primaryNoHit");

		sb.append('\t').append("WKID");
		sb.append('\t').append("KID");
		sb.append('\t').append("ANI");
		sb.append('\t').append("Complt");
		sb.append('\t').append("Contam");
		sb.append('\t').append("TaxID");
		sb.append('\t').append("TaxName");
		sb.append('\t').append("topContamID");
		sb.append('\t').append("topContamName");
		
		sb.append('\n');
		
		return sb.toString();
	}
	
	/**
	 * Represents summary data for a single query's sketch results.
	 * Contains query metadata and list of result lines with methods
	 * for parsing headers and generating formatted output.
	 */
	private class SketchResultsSummary {
		
		/** Constructs summary from query header line.
		 * @param line Query header line containing metadata */
		SketchResultsSummary(String line){
			parseHeader(line);
		}

		/**
		 * Parses query header line to extract metadata fields.
		 * Processes tab-delimited key:value pairs for query name,
		 * sequence count, bases, genome size, and sketch length.
		 * @param line Header line to parse
		 */
		void parseHeader(String line){
			String[] split=line.split("\t");
			for(String s : split){
				String[] split2=s.trim().split(": ");
				assert(split2.length==2) : "\n"+line+"\n"+s+"\n"+Arrays.toString(split2)+"\n";
				String a=split2[0], b=split2[1];
//				outstream.println(a+", "+b);
				if(a.equals("Query")){
					query=b;
				}else if(a.equals("Seqs")){
					seqs=Integer.parseInt(b);
				}else if(a.equals("Bases")){
					bases=Long.parseLong(b);
				}else if(a.equals("gSize")){
					gSize=Long.parseLong(b);
				}else if(a.equals("SketchLen")){
					sketchLen=Integer.parseInt(b);
				}else if(a.equals("TaxID")){
					taxID=Integer.parseInt(b);
				}else if(a.equals("IMG")){
					img=Long.parseLong(b);
				}else if(a.equals("File")){
					sketchLen=Integer.parseInt(b);
				}
			}
		}
		
		/** Adds a result line to this summary.
		 * @param line Result line containing match statistics */
		public void add(String line) {
			SketchResultsLine srl=new SketchResultsLine(line);
			list.add(srl);
		}
		
		@Override
		public String toString(){
			StringBuilder sb=new StringBuilder();
			
			sb.append(query);

			sb.append('\t').append(seqs);
			sb.append('\t').append(bases);
			sb.append('\t').append(gSize);
			sb.append('\t').append(sketchLen);
			
			int primaryHits=0;
			int primaryUnique=0;
			int primaryNoHit=0;

			float WKID=0;
			float KID=0;
			float ANI=0;
			float Complt=0;
			float Contam=0;
			int TaxID=0;
			String TaxName=".";
			int topContamID=0;
			String topContamName=".";
			
			SketchResultsLine first=list.size()>0 ? list.get(0) : null;
			SketchResultsLine second=list.size()>1 ? list.get(1) : null;
			for(int i=2; tree!=null && i<list.size() && failsLevelFilter(first.taxID, second.taxID); i++){
				second=list.get(i);
			}
			if(second!=null && failsLevelFilter(first.taxID, second.taxID)){second=list.get(1);}
			
			if(second!=null && uniqueHitsForSecond){
				for(int i=1; i<list.size(); i++){
					
					SketchResultsLine line=list.get(i);
					if(!failsLevelFilter(first.taxID, line.taxID) && line.unique>second.unique && line.unique>=minUniqueHits){
						second=line;
					}
				}
			}
			
			if(first!=null){
				primaryHits=first.matches;
				primaryUnique=first.unique;
				primaryNoHit=first.noHit;

				WKID=first.wkid;
				KID=first.kid;
				ANI=first.ani;
				Complt=first.complt;
				Contam=first.contam;
				TaxID=first.taxID;
				TaxName=first.name;
			}
			if(second!=null){
				topContamID=second.taxID;
				topContamName=second.name;
			}
			
			sb.append('\t').append(primaryHits);
			sb.append('\t').append(primaryUnique);
			sb.append('\t').append(primaryNoHit);

			sb.append('\t').append(Tools.format("%.2f", WKID));
			sb.append('\t').append(Tools.format("%.2f", KID));
			sb.append('\t').append(Tools.format("%.2f", ANI));
			sb.append('\t').append(Tools.format("%.2f", Complt));
			sb.append('\t').append(Tools.format("%.2f", Contam));
			sb.append('\t').append(TaxID);
			sb.append('\t').append(TaxName);
			sb.append('\t').append(topContamID);
			sb.append('\t').append(topContamName);
			
			sb.append('\n');
			
			return sb.toString();
		}
		
		/**
		 * Tests if two taxonomic IDs fail the level filter.
		 * Finds common ancestor and checks if its taxonomic level
		 * is at or above the configured filter level.
		 * @param a First taxonomic ID
		 * @param b Second taxonomic ID
		 * @return True if the taxa are too closely related
		 */
		private boolean failsLevelFilter(int a, int b) {
			if(a<1 || b<1 || tree==null){return false;}
			int c=tree.commonAncestor(a, b);
			TaxNode tn=tree.getNode(c);
			while(!tn.cellularOrganisms() && tn.levelExtended==TaxTree.NO_RANK_E){tn=tree.getNode(tn.pid);}
			
			return tn.levelExtended<=taxLevel;
		}

		/** Query sequence name or file path */
		String query;
		/** Input file name */
		String fname;
		/** Number of sequences in query */
		int seqs;
		/** Total bases in query sequences */
		long bases;
		/** Estimated genome size */
		long gSize;
		/** Length of sketch (number of hashes) */
		int sketchLen;
		/** NCBI taxonomic identifier */
		int taxID;
		/** IMG identifier */
		long img;
		
		/** List of result lines for this query */
		ArrayList<SketchResultsLine> list=new ArrayList<SketchResultsLine>();
		
	}
	
	/**
	 * Represents a single result line from sketch comparison.
	 * Contains match statistics, taxonomic information, and
	 * similarity metrics for one database match.
	 */
	private class SketchResultsLine{
		
		/**
		 * Parses result line into component fields.
		 * Strips ANSI color codes and parses tab-delimited values
		 * for similarity metrics and match statistics.
		 * @param line Tab-delimited result line
		 */
		SketchResultsLine(String line){
			//Handle colors
			if(line.startsWith(Colors.esc)){
				int first=line.indexOf('m');
				int last=line.lastIndexOf(Colors.esc);
				line=line.substring(first+1, last);
			}
			String[] split=line.replaceAll("%", "").split("\t");
			wkid=Float.parseFloat(split[0]);
			kid=Float.parseFloat(split[1]);
			ani=Float.parseFloat(split[2]);
			complt=Float.parseFloat(split[3]);
			contam=Float.parseFloat(split[4]);
			
			matches=Integer.parseInt(split[5]);
			unique=Integer.parseInt(split[6]);
			noHit=Integer.parseInt(split[7]);
			taxID=Integer.parseInt(split[8]);
			gSize=Integer.parseInt(split[9]);
			gSeqs=Integer.parseInt(split[10]);
			
			name=split[11];
			if(name.equals(".") && split.length>11){
				name=split[12];
			}
		}
		
		/** Weighted k-mer identity percentage */
		float wkid;
		/** K-mer identity percentage */
		float kid;
		/** Average nucleotide identity percentage */
		float ani;
		/** Completeness percentage */
		float complt;
		/** Contamination percentage */
		float contam;
		/** Number of matching k-mers */
		int matches;
		/** Number of unique matching k-mers */
		int unique;
		/** Number of k-mers with no database hit */
		int noHit;
		/** NCBI taxonomic identifier for match */
		int taxID;
		/** Genome size of database match */
		int gSize;
		/** Number of sequences in database match genome */
		int gSeqs;
		/** Taxonomic name of database match */
		String name;
	}
	
	/** Input file paths */
	final ArrayList<String> in;
	/** Output file path */
	final String out;
	
	/** Taxonomic tree for filtering related taxa */
	TaxTree tree=null;
	/** Minimum taxonomic level for contamination detection */
	int taxLevel=TaxTree.GENUS_E;
	/** Whether to select second hit based on unique k-mers */
	boolean uniqueHitsForSecond=false;
	/** Minimum unique hits required for second match consideration */
	int minUniqueHits=3;
	/** Whether to include header line in output */
	boolean printHeader=true;
	
	/** Legacy code from SealStats */
	boolean ignoreSameTaxa=false;
	/** Legacy parameter for ignoring same barcode matches */
	boolean ignoreSameBarcode=false;
	/** Legacy parameter for ignoring same location matches */
	boolean ignoreSameLocation=false;
	/** Legacy parameter for using total as denominator */
	boolean totalDenominator=false;
	/** Legacy parameter for printing total statistics */
	boolean printTotal=true;
	
	/** Output stream for messages */
	PrintStream outstream=System.err;
	
}
