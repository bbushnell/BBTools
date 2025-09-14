package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import gff.GffLine;
import shared.LineParser4;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Trims contigs to eliminate uncovered areas.
 * @author Brian Bushnell
 * @date September 28, 2024
 *
 */
public class TrimContigs {

	public static void main(String[] args){
		Timer t=new Timer();
		TrimContigs x=new TrimContigs(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public TrimContigs(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		Shared.capBufferLen(20);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("coverage") || a.equals("cov") || a.equals("ranges") || a.equals("covranges")){
				covRanges=b;
			}else if(a.equals("minc") || a.equals("mincov") || a.equals("mincoverage")){
				minCoverage=Double.parseDouble(b);
			}else if(a.equals("minp") || a.equals("minpercent")){
				minCoveredPercent=Double.parseDouble(b);
			}else if(a.equals("minl") || a.equals("minlen") || a.equals("minlength")){
				minLength=Integer.parseInt(b);
			}else if(a.equals("trim") || a.equals("trimends") || a.equals("trimmin") || a.equals("mintrim")){
				trimEndsMin=Tools.max(Integer.parseInt(b), 0);
			}else if(a.equals("trimmax") || a.equals("maxtrim")){
				if(Tools.startsWithLetter(b)) {
					if(b.equalsIgnoreCase("big") || b.equalsIgnoreCase("2b") || b.equalsIgnoreCase("2g")) {
						trimEndsMax=2000000000;
					}else {
						throw new RuntimeException("Unknown argument "+arg);
					}
				}else {
					trimEndsMax=(int)Tools.min(Long.parseLong(b), 2000000000);
				}
			}else if(a.equals("trimextra") || a.equals("extra")){
				trimEndsExtra=Integer.parseInt(b);
			}else if(a.equals("maxuncovered")){
				maxUncoveredLength=Integer.parseInt(b);
			}else if(a.equals("break") || a.equals("breakcontigs")){
				breakContigs=Parse.parseBoolean(b);
			}else if(a.equals("breaklist")){
				breakListFile=b;
			}else if(a.equalsIgnoreCase("skippolyn")){
				skipPolyN=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("gff")){
				throw new RuntimeException("Please specify gffin or gffout.");
			}else if(a.equalsIgnoreCase("gffin")){
				gffIn=b;
			}else if(a.equalsIgnoreCase("gffout")){
				gffOut=b;
			}
			
			else if(a.equals("appendresults") || a.equals("logappend") || a.equals("appendlog") || a.equals("appendtolog")){
				logappend=Parse.parseBoolean(b);
			}else if(a.equals("log") || a.equals("results")){
				logfile=b;
			}else if(a.equals("logheader")){
				logheader=Parse.parseBoolean(b);
			}else if(a.equals("outd") || a.equals("outdirty") || a.equals("outb") || a.equals("outbad")){
				outdirty=b;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			if(parser.minReadLength>0){minLength=parser.minReadLength;}
			
			in1=parser.in1;
			qfin1=parser.qfin1;

			outclean=parser.out1;
			qfoutclean=parser.qfout1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		minLength=Tools.max(1, minLength);
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		name=ReadWrite.stripToCore(in1);
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(outclean!=null && outclean.equalsIgnoreCase("null")){outclean=null;}
		if(outdirty!=null && outdirty.equalsIgnoreCase("null")){outdirty=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, outclean, outdirty)){
			outstream.println((outclean==null)+", "+outclean+", "+(outdirty==null)+", "+outdirty);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outclean+", "+outdirty+"\n");
		}

		ffoutclean=FileFormat.testOutput(outclean, FileFormat.FASTA, extout, true, overwrite, append, false);
		ffoutdirty=FileFormat.testOutput(outdirty, FileFormat.FASTA, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		ffRange=FileFormat.testInput(covRanges, FileFormat.TEXT, ".txt", true, false);
		
		assert(ffRange!=null) : "No coverage file specified.";
	}
	
	void process(Timer t){

		final HashMap<String, ArrayList<Range>> map=new HashMap<String, ArrayList<Range>>(1024);
		final LineParser4 lp=new LineParser4("-\t");
		if(ffRange!=null){
			ByteFile tf=ByteFile.makeByteFile(ffRange);
			int i=0;
			String name=null;
			ArrayList<Range> list=null;
			for(byte[] s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s[0]=='#'){
					name=new String(s, 1, s.length-1);
					list=new ArrayList<Range>(2);
					Object old=map.put(name, list);
					assert(old==null) : "Duplicate contig name "+name;
				}else{
					lp.set(s);
					int a=lp.parseInt(0);
					int b=lp.parseInt(1);
					float d=lp.parseFloat(2);
					list.add(new Range(a, b, d));
				}
				i++;
			}
			tf.close();
		}
		
		if(gffIn!=null) {loadGff(gffIn);}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, qfin1, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		assert(!cris.paired());
		
		final ConcurrentReadOutputStream rosClean;
		if(outclean!=null){
			final int buff=4;
			
			assert(!outclean.equalsIgnoreCase(in1) && !outclean.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			rosClean=ConcurrentReadOutputStream.getStream(ffoutclean, null, qfoutclean, null, buff, null, false);
			rosClean.start();
		}else{rosClean=null;}
		
		final ConcurrentReadOutputStream rosDirty;
		if(outdirty!=null){
			final int buff=4;
			
			assert(!outdirty.equalsIgnoreCase(in1) && !outdirty.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			rosDirty=ConcurrentReadOutputStream.getStream(ffoutdirty, null, qfoutdirty, null, buff, null, false);
			rosDirty.start();
		}else{rosDirty=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;

		long basesTrimmed=0;
		long contigsBroken=0;
		long contigBreaks=0;
		long contigsTrimmed=0;
		
		long readsOut=0;
		long basesOut=0;
		
		long readsFiltered=0;
		long basesFiltered=0;
		
		long inputSeqs=0, outputSeqs=0, dirtySeqs=0;
		long inputBases=0, outputBases=0;
		long dirtyBases=0;
		
		ByteStreamWriter bswBreaks=(breakListFile==null ? null : 
			ByteStreamWriter.makeBSW(breakListFile, overwrite, append, true));
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				final ArrayList<Read> cleanList=new ArrayList<Read>(reads.size());
				final ArrayList<Read> dirtyList=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read seq=reads.get(idx);
					assert(seq.mate==null);
					
					final int initialLength1=seq.length();
					
					readsProcessed++;
					basesProcessed+=initialLength1;
					
					final ArrayList<Range> ranges=map.get(seq.id);
					if(ranges!=null && ranges.size()>1) {
						if(!breakContigs) {
							Range r=toMaximalRange(ranges);
							ranges.clear();
							ranges.add(r);
						}else if(skipPolyN || maxUncoveredLength>0) {
							fixPolyN(seq, ranges);
						}
					}
					ArrayList<GffLine> gffLines=(gffMap==null ? null : gffMap.get(seq.id));
					
					if(ranges==null || ranges.isEmpty() || seq.length()-2*trimEndsMin<minLength){//Common case
//						System.err.println("Discarding bad contig "+seq.name());
						dirtyList.add(seq);
					}else if(ranges.size()==1){//Common case
						Range r=ranges.get(0);
						r.name=seq.id;
						Read processed=processSeq(seq, r, 1, ranges.size());
						if(processed.discarded()) {dirtyList.add(processed);}
						else {
							contigsTrimmed+=processed.length()<seq.length() ? 1 : 0;
							cleanList.add(processed);
							processGff(ranges, gffLines, gffLinesOut);
						}
						contigsTrimmed+=(processed.length()==seq.length() ? 0 : 1);
					}else {//Uncommon
						contigsTrimmed++;
						contigsBroken++;
						for(int i=0; i<ranges.size(); i++) {
							Range r=ranges.get(i);
							Read processed=processSeq(seq, r, i, ranges.size());
							contigBreaks++;
							if(processed==null) {
								//This doesn't actually happen
							}else {
								r.name=processed.name();
								if(processed.discarded()) {dirtyList.add(processed);}
								else {cleanList.add(processed);}
							}
						}
						processGff(ranges, gffLines, gffLinesOut);
						contigBreaks--;
//						contigBreakList.add(seq.id);
						if(bswBreaks!=null) {bswBreaks.println(seq.id);}
					}
				}
				
				for(Read r : reads) {
					inputSeqs++;
					inputBases+=r.length();
				}
				for(Read r : cleanList) {
					outputSeqs++;
					outputBases+=r.length();
				}
				for(Read r : dirtyList) {
					dirtySeqs++;
					dirtyBases+=r.length();
				}
				
				if(rosClean!=null){rosClean.add(cleanList, ln.id);}
				if(rosDirty!=null){rosDirty.add(dirtyList, ln.id);}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		readsOut=outputSeqs;
		basesOut=outputBases;
		readsFiltered=dirtySeqs;
		basesFiltered=dirtyBases;
		basesTrimmed=(inputBases-(basesOut+basesFiltered));
		
		if(gffOut!=null) {
			ByteStreamWriter bswGff=ByteStreamWriter.makeBSW(gffOut, overwrite, append, true);
			bswGff.println("##gff-version 3");
			final ByteBuilder bb=new ByteBuilder();
			for(GffLine line : gffLinesOut) {
				line.appendTo(bb).nl();
				bswGff.print(bb);
				bb.clear();
			}
			bswGff.poisonAndWait();
		}
		errorState=ReadWrite.closeStreams(cris, rosClean, rosDirty)||errorState;
		if(bswBreaks!=null) {errorState=bswBreaks.poisonAndWait()||errorState;}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);
		
		outstream.println("Time:               "+t);
		outstream.println("Scaffolds In:       "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases In:           "+basesProcessed+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Scaffolds Out:      "+readsOut);
		outstream.println("Bases Out:          "+basesOut);
		outstream.println("Scaffolds Filtered: "+dirtySeqs);
		outstream.println("Bases Filtered:     "+dirtyBases);
		outstream.println("Scaffolds Trimmed:  "+contigsTrimmed);
		outstream.println("Bases Trimmed:      "+basesTrimmed);
		outstream.println("Scaffolds Broken:   "+contigsBroken);
		outstream.println("Scaffold Breaks:    "+contigBreaks);
//		outstream.println("Broken Scaffolds:   "+contigBreakList);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	private Read processSeq(Read seq, Range r, int partnum, int parts) {
//		float cov=(r.depth*(r.b-r.a+1))/(seq.length());
		float cov=r.depth;
		r.name=seq.id;
		if(parts==1) {//Common case
			if(cov>=minCoverage && seq.length()>=minLength && trimEndsMin<1 
					&& r.a<=maxUncoveredLength && seq.length()-r.b-1<=maxUncoveredLength) {
				//Hopefully, common case.
//				System.err.println("1) Returning good contig "+seq.name());
				return seq;
			}else if(seq.length()<minLength || cov<minCoverage) {
				//Another common case
				seq.setDiscarded(true);
//				System.err.println("2) Returning bad contig "+seq.name());
				return seq;
			}else if(r.length()+2*maxUncoveredLength<minLength) {
//				System.err.println("3) Returning bad contig "+seq.name());
				//Another common case
				seq.setDiscarded(true);
				return seq;
			}
		}
		
		if(r.a>=maxUncoveredLength) {r.a+=trimEndsExtra;}
		else {r.a=0;}
		r.a=Tools.mid(r.a, trimEndsMin, trimEndsMax);
		if(seq.length()-r.b-1>maxUncoveredLength) {r.b-=trimEndsExtra;}
		else {r.b=seq.length()-1;}
		r.b=Tools.mid(r.b, seq.length()-trimEndsMin-1, seq.length()-trimEndsMax-1);
		
		int trimmed=0;
		if(r.a>0 || r.b<seq.length()-1) {
			seq=seq.clone();
			trimmed=TrimRead.trimToPosition(seq, r.a, r.b, 0);
			if(parts>1) {seq.id=seq.id+"_part"+partnum;}
		}
		
		seq.setDiscarded(cov<minCoverage || seq.length()<minLength);
		r.setDiscarded(seq.discarded());
//		System.err.println("4) Returning "+(seq.discarded() ? "bad" : "good")+
//				" contig "+seq.name()+", trimmed="+trimmed);
		r.name=seq.id;
		return seq;
	}
	
	private static void processGff(ArrayList<Range> ranges, 
			ArrayList<GffLine> gffLines, ArrayList<GffLine> gffLinesOut) {
		if(ranges==null || gffLines==null || gffLines.isEmpty()) {return;}
		for(Range r :ranges) {
//			System.err.println("Processing "+r);
			if(!r.discarded) {
				int rlen=r.length();
				for(GffLine line : gffLines) {
//					System.err.println("Processing "+line);
					int start=Tools.max(r.a+1, line.start);
					int stop=Tools.min(r.b+1, line.stop);
					if(start<=stop) {//Overlap
						int a=start-1, b=stop-1;
						int a2=a-r.a;
						int b2=b-r.a;
						int start2=a2+1, stop2=b2+1;
						GffLine line2=line.clone();
						line2.seqid=r.name;
						line2.start=start2;
						line2.stop=stop2;
						
						if(line2.strand()==GffLine.PLUS) {
							if(a<r.a && line2.phase>=0) {
								int trimmed=r.a-a;
								line2.phase=(line2.phase+trimmed)%3;
							}
						}else {
							if(b>r.b && line2.phase>=0) {
								int trimmed=b-r.b;
								line2.phase=(line2.phase+trimmed)%3;
							}
						}
						
						if(line2.length()>=line.length()/4) {gffLinesOut.add(line2);}
					}
				}
			}
		}
	}
	
	private static Range toMaximalRange(Collection<Range> ranges) {
		assert(ranges.size()>1) : ranges.size();
		int min=Integer.MAX_VALUE;
		int max=-1;
		double depthSum=0;
		for(Range r : ranges) {
			depthSum+=r.depthSum();
			min=Tools.min(min, r.a);
			max=Tools.max(max, r.b);
		}
		Range r=new Range(min, max);
		r.depth=(float)(depthSum/r.length());
		return r;
	}
	
	/** Fuses adjacent ranges separated by fewer than 
	 * maxUncoveredLength uncovered defined bases. */
	private int fixPolyN(Read seq, ArrayList<Range> ranges) {
		assert(ranges.size()>1);
		final byte[] bases=seq.bases;
		ArrayList<Range> fixed=new ArrayList<Range>(ranges.size());
		int fusions=0;
		for(int i=0, j=1; j<ranges.size(); i++, j++) {
			Range left=ranges.get(i), right=ranges.get(j);
			int defined=0, undefined=0;
			for(int k=left.b+1; k<right.a; k++) {
				byte b=bases[k];
				if(AminoAcid.isFullyDefined(b)) {
					defined++;
				}else {
					undefined++;
				}
			}
			if(!skipPolyN) {
				defined+=undefined;
				undefined=0;
			}
			if(defined<=maxUncoveredLength || 
					(undefined>0 && defined<=maxUncoveredLength*2)) {//fuse
				right.absorb(left);
				fusions++;
			}else {
				fixed.add(left);
			}
		}
		fixed.add(ranges.get(ranges.size()-1));
		ranges.clear();
		ranges.addAll(fixed);
		return fusions;
	}
	
	/*--------------------------------------------------------------*/
	
	private void loadGff(String s) {
		gffLinesIn=GffLine.loadGffFile(gffIn, null, false);
		gffMap=new HashMap<String, ArrayList<GffLine>>();
		for(GffLine line : gffLinesIn) {
			ArrayList<GffLine> list=gffMap.get(line.seqid);
			if(list==null) {
				list=new ArrayList<GffLine>();
				gffMap.put(line.seqid, list);
			}
			list.add(line);
		}
		gffLinesOut=new ArrayList<GffLine>();
//		assert(false) : gffMap.size()+", "+gffLinesIn.size();
	}
	
	private static class Range {
		public Range(int a_, int b_) {
			a=a_;
			b=b_;
		}
		public Range(int a_, int b_, float depth_) {
			a=a_;
			b=b_;
			depth=depth_;
		}
		public void absorb(Range r) {
			float depthsum=depthSum()+r.depthSum();
			a=Tools.min(a, r.a);
			b=Tools.max(b, r.b);
			depth=(depth==-1 ? -1 : depthsum/length());
		}
		public void setDiscarded(boolean discarded_) {
			discarded=discarded_;
		}
		@Override
		public String toString() {
			return name+"("+a+"-"+b+")";
		}
		public int length() {return b-a+1;}
		public float depthSum() {return depth*length();}
		int a;
		int b;
		float depth=-1;
		boolean discarded=false;
		String name;
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String covRanges=null;
	private String name=null;
	
	private String qfin1=null;

	private String outclean=null;
	private String outdirty=null;

	private String qfoutclean=null;
	private String qfoutdirty=null;
	
	private String extin=null;
	private String extout=null;
	
	private String gffIn=null;
	private String gffOut=null;
	
	/*--------------------------------------------------------------*/

	private ArrayList<GffLine> gffLinesIn;
//	private HashMap<StringNum, GffLine> lineMap;
	private HashMap<String, ArrayList<GffLine>> gffMap;
	private ArrayList<GffLine> gffLinesOut;
	
	/*--------------------------------------------------------------*/
	
	private long maxReads=-1;

	/** Scaffolds shorter than this will be discarded. */
	private int minLength=1;
	/** Scaffolds with lower average coverage will be discarded. */
	private double minCoverage=1;
	/** Scaffolds with a lower percent of covered bases will be discarded. */
	private double minCoveredPercent=0;

	/** Trim this much from sequence ends, minimum */
	private int trimEndsMin=0;
	
	/** Trim this much from sequence ends, maximum */
	private int trimEndsMax=2000000000;
	
	/** Trim this far into covered areas */
	private int trimEndsExtra=5;
	
	/** Uncovered areas longer than this will trigger trimming */
	private int maxUncoveredLength=3;
	
	/** Permission to break apart contigs at uncovered areas */
	boolean breakContigs=true;
	boolean skipPolyN=true;
	String breakListFile;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
//	private final FileFormat ffCov;
	private final FileFormat ffRange;

	private final FileFormat ffoutclean;
	private final FileFormat ffoutdirty;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean logappend=false;
	private String logfile=null;
	private boolean logheader=true;
	private static boolean PRINT_SHORT_CONTIG_RESULTS=false;
	
}
