package barcode;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import parse.LineParserS1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.LongList;

/**
 * Compares the last two delimiter-separated fields of each primary read's header.
 * Counts one label pair per fragment, while the read-throughput counter includes
 * mates. UNKNOWN is recognized case-insensitively; other label equality and quantSet
 * membership are case-sensitive. Mutable parsing/counter state has one processing owner.
 * Normalized reports require observed labels; empty runs without reports are allowed.
 * @author Brian Bushnell
 * @date May 8, 2014
 */
public class CompareLabels{

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CompareLabels x=new CompareLabels(args);
		
		//Close redirected diagnostics even when processing reports a stream failure.
		try{x.process(t);}finally{Shared.closeStream(x.outstream);}
	}
	
	/**
	 * Constructs a CompareLabels instance with command-line arguments.
	 * Parses input/output files, delimiter settings, and filtering options.
	 * @param args Command-line arguments including input files and parameters
	 */
	public CompareLabels(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		parser.out1=out1;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("delimiter")){
				delimiter=Parse.parseSymbolToCharacter(b);
			}else if(a.equals("swap")){
				swap=Parse.parseBoolean(b);
			}else if(a.equals("labelstats")){
				labelStats=b;
			}else if(a.equals("quantset")){
				quantSetFile=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		lp=new LineParserS1(delimiter);
		//BAR213: apply parsed output policy to both report writers.
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		
		map=(labelStats==null ? null : new LinkedHashMap<String, Label>());
		if(quantSetFile!=null){
			quantSet=new HashSet<String>();
			quantSet.add("UNKNOWN");
			for(String s : TextFile.toStringLines(quantSetFile)){quantSet.add(s);}
		}else{
			quantSet=null;
		}
	}
	
	/**
	 * Main processing method that compares labels and generates statistics.
	 * Reads input sequences with dual labels, computes agreement/disagreement rates, and outputs summary statistics and per-label breakdown.
	 * @param t Timer for tracking execution time
	 * @throws RuntimeException If input or output shutdown reports an error
	 * @throws IllegalStateException If a normalized report is requested without observed labels
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}
		boolean paired=cris.paired();
		
		processInner(cris);
		
		//BAR212: retain failure flags while completing both report writers.
		boolean errorState=ReadWrite.closeStream(cris);
		//BAR218: reject undefined normalized reports before either output is opened.
		if(labelsProcessed==0 && (ffout1!=null || labelStats!=null)){
			throw new IllegalStateException("CompareLabels cannot write normalized reports without observed labels.");
		}
		if(verbose){outstream.println("Finished.");}
		
		if(ffout1!=null){
			ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);
			summarize(bsw);
			errorState|=bsw.poisonAndWait();
		}
		if(labelStats!=null){
			FileFormat ff=FileFormat.testOutput(labelStats, FileFormat.TXT, null, true, overwrite, append, false);
			ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
			printMap(bsw);
			errorState|=bsw.poisonAndWait();
		}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		if(errorState){throw new RuntimeException("CompareLabels encountered an input or output stream error; output may be incomplete.");}
	}
	
	/**
	 * Inner processing loop that handles batched read processing.
	 * Fetches read batches from input stream and processes each read's labels.
	 * @param cris Concurrent read input stream for efficient parallel processing
	 */
	private void processInner(ConcurrentReadInputStream cris){
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
		}

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				readsProcessed+=r1.pairCount();
				basesProcessed+=r1.pairLength();
				labelsProcessed++;

				processRead(r1);
			}

			cris.returnList(ln);
			if(verbose){outstream.println("Returned a list.");}
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	/**
	 * Generates summary statistics comparing label agreement rates.
	 * Computes yields and contamination fractions from the five label-pair categories.
	 * Label totals include invalid pairs; assigned totals exclude UNKNOWN and invalid pairs.
	 * Contamination fields and per-category PPM columns use parts per million.
	 * @param bsw Output writer for summary statistics (null to skip output)
	 */
	private void summarize(ByteStreamWriter bsw){
		if(bsw==null){return;}
		ByteBuilder bb=new ByteBuilder();
		long valid=labelsProcessed-invalidCount;
		float fracMult=1f/(Tools.max(1, labelsProcessed));
		float ppmMult=1000000f/labelsProcessed;
		

		long count1=aaCount+auCount+abCount;
		long count2=aaCount+uaCount+abCount;
		float mult1=1f/Tools.max(count1, 1);
		float mult2=1f/Tools.max(count2, 1);
		float ryield1=aaCount*mult2;
		float ryield2=aaCount*mult1;
		float ayield1=count1*fracMult;
		float ayield2=count2*fracMult;
		float contam1=abCount*mult1;
		float contam2=abCount*mult2;
		
		bb.append("#Labels").tab().append(labelsProcessed).nl();
		bb.append("#Valid").tab().append(valid).tab().append(fracMult*valid, 6).nl();
		bb.append("#TermsPerRead").tab().append(termCounts.meanHist(), 6).nl();
		bb.append("#RelYield1").tab().append(ryield1, 5).nl();
		bb.append("#RelYield2").tab().append(ryield2, 5).nl();
		bb.append("#AbsYield1").tab().append(ayield1, 5).nl();
		bb.append("#AbsYield2").tab().append(ayield2, 5).nl();
		//BAR206: convert contamination fractions to the units named by these fields.
		bb.append("#Contam1_PPM").tab().append(contam1*1000000, 2).nl();
		bb.append("#Contam2_PPM").tab().append(contam2*1000000, 2).nl();
		bb.append("#Metric").tab().append("Count").tab().append("Rate").tab().append("PPM").nl();
		bb.append("AACount").tab().append(aaCount).tab().append(fracMult*aaCount, 5).tab().append(ppmMult*aaCount, 2).nl();
		bb.append("UUCount").tab().append(uuCount).tab().append(fracMult*uuCount, 5).tab().append(ppmMult*uuCount, 2).nl();
		bb.append("AUCount").tab().append(auCount).tab().append(fracMult*auCount, 5).tab().append(ppmMult*auCount, 2).nl();
		bb.append("UACount").tab().append(uaCount).tab().append(fracMult*uaCount, 5).tab().append(ppmMult*uaCount, 2).nl();
		bb.append("ABCount").tab().append(abCount).tab().append(fracMult*abCount, 5).tab().append(ppmMult*abCount, 2).nl();
		bsw.println(bb);
	}
	
	/**
	 * Outputs per-label statistics sorted by occurrence frequency.
	 * Shows individual label performance including agreement counts, yields, and contamination.
	 * Places UNKNOWN first, then orders each group by descending count and lexical name.
	 * @param bsw Output writer for per-label statistics (null to skip output)
	 */
	private void printMap(ByteStreamWriter bsw){
		if(bsw==null){return;}
		ByteBuilder bb=new ByteBuilder();
		long valid=labelsProcessed-invalidCount;
		float fracMult=1f/labelsProcessed;
		float ppmMult=1000000f/labelsProcessed;
		bb.append("#Labels").tab().append(labelsProcessed).nl();
		bb.append("#Valid").tab().append(valid).tab().append(fracMult*valid, 5).nl();
		bb.append("#TermsPerRead").tab().append(termCounts.meanHist(), 5).nl();
		bb.append(Label.header());
		bsw.println(bb);
		
		ArrayList<Label> list=new ArrayList<Label>();
		for(Entry<String, Label> e : map.entrySet()){
			Label lab=e.getValue();
			list.add(lab);
		}
		Collections.sort(list);
		for(Label lab : list){
			bb.clear();
			lab.appendTo(bb);
			bsw.println(bb);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Processes a single read's dual labels for comparison analysis.
	 * Extracts two labels from read ID, categorizes agreement/disagreement, and updates appropriate counters for statistical analysis.
	 * @param r Read containing dual labels in the sequence identifier
	 */
	private void processRead(Read r){
		lp.set(r.id);
		final int terms=lp.terms();
		termCounts.increment(terms);
		if(terms<2){
			invalidCount++;
			return;
		}
		final String s1=lp.parseString(terms-(swap ? 1 : 2));
		final String s2=lp.parseString(terms-(swap ? 2 : 1));
		final boolean unknown1=s1.equalsIgnoreCase("UNKNOWN");
		final boolean unknown2=s2.equalsIgnoreCase("UNKNOWN");
		final boolean equal=s1.equals(s2);
		if(quantSet!=null && !(quantSet.contains(s1) && quantSet.contains(s2))){
			invalidCount++;
			return;
		}
		if(unknown1){
			if(unknown2){uuCount++;}
			else {uaCount++;}
		}else if(unknown2){
			auCount++;
		}else if(s1.equals(s2)){
			aaCount++;
		}else{
			abCount++;
		}
		
		if(map!=null){
			final Label lab1=getLabel(s1);
			final Label lab2=(equal ? lab1 : getLabel(s2));
			if(equal){lab1.aa++;}
			else if(unknown1){
				lab1.ua++;
				lab2.ua++;
			}else if(unknown2){
				lab1.au++;
				lab2.au++;
			}else{
				lab1.ab++;
				lab2.ba++;
			}
		}
	}
	
	/**
	 * Retrieves or creates a Label object for the given label string.
	 * Maintains a map of all observed labels for per-label statistics tracking.
	 * @param s Label string to look up or create
	 * @return Label object for tracking statistics, or null if input is null
	 */
	private Label getLabel(String s){
		if(s==null){return null;}
		Label l=map.get(s);
		if(l==null){
			l=new Label(s);
			map.put(s, l);
		}
		return l;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Per-name counts: AA agreement, AU/UA one unknown, AB/BA differing assigned labels in each direction. */
	private static class Label implements Comparable<Label>{
		
		/** Retains the label name and caches case-insensitive UNKNOWN classification. */
		Label(String s){
			name=s;
			unknown="UNKNOWN".equalsIgnoreCase(name);
		}
		
		/** Total contributions to this label's row; a disagreement can contribute to two distinct rows. */
		final long count(){return aa+au+ua+ab+ba;}
		
		//aka refcount
		final long count1(){return aa+ab+(unknown ? ua : au);}
		
		//aka assigned count
		final long count2(){return aa+ba+(unknown ? au : ua);}

		final long tp(){return aa;}
		final long fp(){return ba;}
		final long fn(){return au;}
		
		/** Orders UNKNOWN before ordinary labels, then descending count and lexical name within each group. */
		@Override
		public int compareTo(Label o){
			//BAR205 / barcode/CompareLabels#001: apply UNKNOWN priority symmetrically in both argument orders.
			if(unknown!=o.unknown){return unknown ? -1 : 1;}
			long dif=count()-o.count();
			if(dif>0){return -1;}
			else if(dif<0){return 1;}
			return name.compareTo(o.name);
		}
		
		/** Column names for the row emitted by appendTo. */
		public static String header(){
			ByteBuilder bb=new ByteBuilder();
			bb.append("#Name").tab().append("Count").tab().append("Count1").tab().append("Count2");
			bb.tab().append("AA").tab().append("AU").tab().append("UA");
			bb.tab().append("AB").tab().append("BA");
			bb.tab().append("Yield1").tab().append("Yield2");
			bb.tab().append("Contam1").tab().append("Contam2");
			return bb.toString();
		}
		
		/**
		 * Appends this label's statistics to a ByteBuilder in tabular format.
		 * Includes counts, fractional yields and contamination in parts per million.
		 * @param bb ByteBuilder to append statistics to
		 * @return The same ByteBuilder for method chaining
		 */
		public ByteBuilder appendTo(ByteBuilder bb){
			long count=count();
			long count1=count1();
			long count2=count2();
			float mult1=1f/Tools.max(count1, 1);
			float mult2=1f/Tools.max(count2, 1);
			float yield1=aa*mult2;
			float yield2=aa*mult1;
			float contam1=ab*mult1;
			float contam2=ba*mult2;
			bb.append(name);
			bb.tab().append(count).tab().append(count1).tab().append(count2);
			bb.tab().append(aa).tab().append(au).tab().append(ua);
			bb.tab().append(ab).tab().append(ba);
			bb.tab().append(yield1, 5).tab().append(yield2, 5);
			bb.tab().append(contam1*1000000, 2).tab().append(contam2*1000000, 2);
			return bb;
		}
		
		final String name;
		final boolean unknown;
		long aa=0;
		long au=0;
		long ua=0;
		long ab=0;
		long ba=0;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1="stdout.txt";
	/** Parsed output policy shared by the summary and per-label reports. */
	private final boolean overwrite, append;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	private char delimiter='\t';
	private final LineParserS1 lp;

	private boolean swap=false;
	private String labelStats=null;
	private String quantSetFile=null;

	private final LinkedHashMap<String, Label> map;
	private final HashSet<String> quantSet;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/** Individual reads and bases, including mates. */
	private long readsProcessed=0;
	private long basesProcessed=0;
	/** Primary-read label pairs, including pairs later classified invalid. */
	private long labelsProcessed=0;
	
	private LongList termCounts=new LongList();
	private long uuCount=0;
	private long uaCount=0;
	private long auCount=0;
	private long aaCount=0;
	private long abCount=0;
	private long invalidCount=0;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
