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
import shared.LineParserS1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date May 8, 2014
 *
 */
public class CompareLabels {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CompareLabels x=new CompareLabels(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
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
		}
		
		lp=new LineParserS1(delimiter);
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		
		map=(labelStats==null ? null : new LinkedHashMap<String, Label>());
		if(quantSetFile!=null) {
			quantSet=new HashSet<String>();
			quantSet.add("UNKNOWN");
			for(String s : TextFile.toStringLines(quantSetFile)) {quantSet.add(s);}
		}else {
			quantSet=null;
		}
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}
		boolean paired=cris.paired();
		
		processInner(cris);
		
		ReadWrite.closeStream(cris);
		if(verbose){outstream.println("Finished.");}
		
		if(ffout1!=null) {
			ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);
			summarize(bsw);
			bsw.poisonAndWait();
		}
		if(labelStats!=null) {
			FileFormat ff=FileFormat.testOutput(labelStats, FileFormat.TXT, null, true, true, false, false);
			ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
			printMap(bsw);
			bsw.poisonAndWait();
		}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}
	
	private void processInner(ConcurrentReadInputStream cris) {
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
	
	private void summarize(ByteStreamWriter bsw) {
		if(bsw==null) {return;}
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
		bb.append("#Contam1_PPM").tab().append(contam1, 2).nl();
		bb.append("#Contam2_PPM").tab().append(contam2, 2).nl();
		bb.append("#Metric").tab().append("Count").tab().append("Rate").tab().append("PPM").nl();
		bb.append("AACount").tab().append(aaCount).tab().append(fracMult*aaCount, 5).tab().append(ppmMult*aaCount, 2).nl();
		bb.append("UUCount").tab().append(uuCount).tab().append(fracMult*uuCount, 5).tab().append(ppmMult*uuCount, 2).nl();
		bb.append("AUCount").tab().append(auCount).tab().append(fracMult*auCount, 5).tab().append(ppmMult*auCount, 2).nl();
		bb.append("UACount").tab().append(uaCount).tab().append(fracMult*uaCount, 5).tab().append(ppmMult*uaCount, 2).nl();
		bb.append("ABCount").tab().append(abCount).tab().append(fracMult*abCount, 5).tab().append(ppmMult*abCount, 2).nl();
		bsw.println(bb);
	}
	
	private void printMap(ByteStreamWriter bsw) {
		if(bsw==null) {return;}
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
		for(Entry<String, Label> e : map.entrySet()) {
			Label lab=e.getValue();
			list.add(lab);
		}
		Collections.sort(list);
		for(Label lab : list) {
			bb.clear();
			lab.appendTo(bb);
			bsw.println(bb);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void processRead(Read r) {
		lp.set(r.id);
		final int terms=lp.terms();
		termCounts.increment(terms);
		if(terms<2) {
			invalidCount++;
			return;
		}
		final String s1=lp.parseString(terms-(swap ? 1 : 2));
		final String s2=lp.parseString(terms-(swap ? 2 : 1));
		final boolean unknown1=s1.equalsIgnoreCase("UNKNOWN");
		final boolean unknown2=s2.equalsIgnoreCase("UNKNOWN");
		final boolean equal=s1.equals(s2);
		if(quantSet!=null && !(quantSet.contains(s1) && quantSet.contains(s2))) {
			invalidCount++;
			return;
		}
		if(unknown1) {
			if(unknown2) {uuCount++;}
			else {uaCount++;}
		}else if(unknown2) {
			auCount++;
		}else if(s1.equals(s2)) {
			aaCount++;
		}else {
			abCount++;
		}
		
		if(map!=null) {
			final Label lab1=getLabel(s1);
			final Label lab2=(equal ? lab1 : getLabel(s2));
			if(equal) {lab1.aa++;}
			else if(unknown1) {
				lab1.ua++;
				lab2.ua++;
			}else if(unknown2) {
				lab1.au++;
				lab2.au++;
			}else {
				lab1.ab++;
				lab2.ba++;
			}
		}
	}
	
	private Label getLabel(String s) {
		if(s==null) {return null;}
		Label l=map.get(s);
		if(l==null) {
			l=new Label(s);
			map.put(s, l);
		}
		return l;
	}
	
	/*--------------------------------------------------------------*/
	
	private static class Label implements Comparable<Label> {
		
		Label(String s) {
			name=s;
			unknown="UNKNOWN".equalsIgnoreCase(name);
		}
		
		final long count() {return aa+au+ua+ab+ba;}
		
		//aka refcount
		final long count1() {return aa+ab+(unknown ? ua : au);}
		
		//aka assigned count
		final long count2() {return aa+ba+(unknown ? au : ua);}

		final long tp() {return aa;}
		final long fp() {return ba;}
		final long fn() {return au;}
		
		@Override
		public int compareTo(Label o) {
			//Put unknown at the top
			if(name.equalsIgnoreCase("UNKNOWN") && !o.name.equalsIgnoreCase("UNKNOWN")) {return -1;}
			long dif=count()-o.count();
			if(dif>0) {return -1;}
			else if(dif<0) {return 1;}
			return name.compareTo(o.name);
		}
		
		public static String header() {
			ByteBuilder bb=new ByteBuilder();
			bb.append("#Name").tab().append("Count").tab().append("Count1").tab().append("Count2");
			bb.tab().append("AA").tab().append("AU").tab().append("UA");
			bb.tab().append("AB").tab().append("BA");
			bb.tab().append("Yield1").tab().append("Yield2");
			bb.tab().append("Contam1").tab().append("Contam2");
			return bb.toString();
		}
		
		public ByteBuilder appendTo(ByteBuilder bb) {
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
	
	private long readsProcessed=0;
	private long basesProcessed=0;
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
