package jgi;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.NameMapper;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import structures.ListNum;

/**
 * Loads reference sequences from {@code in=} and two primer-alignment maps from
 * SAM files, then emits one subsequence for each reference whose primer spans
 * are both present and non-overlapping.
 * <p>
 * The SAM maps are keyed by {@code RNAME} and looked up by exact input read id
 * or an unambiguous first-token alias. With {@code include=f}, the output omits
 * both primers and keeps only the bases between them. With {@code include=t},
 * the output keeps the clamped primer-to-primer span. Reads without an extracted
 * interval are either dropped or, when {@code fake=t}, represented by a single
 * {@code N} base.
 *
 * @author Brian Bushnell
 * @date Nov 24, 2014
 */
public class CutPrimers {

	/**
	 * Parses arguments, processes configured input, and closes the possibly
	 * redirected status stream.
	 *
	 * @param args command-line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		CutPrimers x=new CutPrimers(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses configuration, generic read options, and local primer options, then
	 * creates input and output format descriptors for {@link #process(Timer)}.
	 *
	 * @param args command-line arguments
	 */
	public CutPrimers(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("sam1")){
				sam1=b;
			}else if(a.equals("sam2")){
				sam2=b;
			}else if(a.equals("fake") || a.equals("addfake")){
				ADD_FAKE_READS=Parse.parseBoolean(b);
			}else if(a.equals("include") || a.equals("includeprimer") || a.equals("includeprimers")){
				INCLUDE_PRIMERS=Parse.parseBoolean(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
	}
	
	/**
	 * Loads the primer SAM maps, streams input reads, writes one output list for
	 * each input list when output is configured, and on loop completion closes
	 * streams and prints input and extraction counts.
	 *
	 * @param t timer to stop and report after processing
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}

		// toSamLines keys each primer alignment by SAM RNAME for input-id lookup.
		NameMapper<SamLine> p1set=toSamLines(sam1);
		NameMapper<SamLine> p2set=toSamLines(sam2);
		long readsProcessed=0, readsSuccess=0;
		RuntimeException processingFailure=null;
		try{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				ArrayList<Read> readsOut=new ArrayList<Read>(reads.size());
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					readsProcessed++;
					final Read r=reads.get(idx);

					// Primer rows are keyed by RNAME; only matching input ids can trim.
					final SamLine sl1=p1set.get(r.id);
					final SamLine sl2=p2set.get(r.id);
					
					// Remember the count to detect whether extraction emits an interval.
					int oldSize=readsOut.size();
					
					final int len=r.length();
					if(sl1!=null && sl2!=null){
						// JGI-129: stop() needs the raw start, including any leading soft clip.
						final int start1=sl1.start(true, false), start2=sl2.start(true, false);
						final int a1=Tools.mid(0, len, start1);
						final int a2=Tools.mid(0, len, start2);
						final int b1=Tools.mid(0, len, sl1.stop(start1, true, false));
						final int b2=Tools.mid(0, len, sl2.stop(start2, true, false));
						if(Tools.overlap(a1, b1, a2, b2)){
							// Overlapping primer spans do not define a retained interval.
						}else{
							
							// Order by clamped start and copy either between or through primers.
							final int from, to;
							if(INCLUDE_PRIMERS){
								// JGI-129: an endpoint at len is already the exclusive copy limit.
								if(a1<a2){
									from=a1;
									to=(b2<len ? b2+1 : len);
								}else{
									from=a2;
									to=(b1<len ? b1+1 : len);
								}
							}else{
								if(a1<a2){
									from=b1+1;
									to=a2;
								}else{
									from=b2+1;
									to=a1;
								}
							}
							
							assert(from>=0 && from<r.bases.length && to>=from) : from+", "+to+", "+r.bases.length+"\n"+
								new String(r.bases)+"\n"+sl1+"\n"+sl2+"\n";
							final byte[] bases=KillSwitch.copyOfRange(r.bases, from, to);
							final byte[] quals=(r.quality==null ? null : KillSwitch.copyOfRange(r.quality, from, to));
							readsOut.add(new Read(bases, quals, r.id, r.numericID));
							readsSuccess++;
						}
					}
					
					// Missing or overlapping primer pairs can still retain output ids.
					if(oldSize==readsOut.size() && ADD_FAKE_READS){
						readsOut.add(new Read(new byte[] {'N'}, null, r.id, r.numericID));
					}
				}
				
				if(ros!=null){ros.add(readsOut, ln.id);}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}catch(RuntimeException e){
			processingFailure=e;
		}
		ReadWrite.closeStreams(cris, ros);
		if(processingFailure!=null){throw processingFailure;}
		if(verbose){outstream.println("Finished.");}
		
		// readsSuccess excludes fake N rows emitted for failed extractions.
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:      "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Sequences Generated:  "+readsSuccess);
	}
	
	/**
	 * Parses non-header SAM rows from one primer file.
	 * <p>
	 * Each retained {@link SamLine} is stored by {@code SamLine.rname()}, replacing
	 * repeated keys according to {@link NameMapper} rules. Input-read lookup later
	 * uses the read id, not the SAM query name.
	 *
	 * @param fname SAM file name
	 * @return primer alignments keyed by mapped reference name
	 */
	public static NameMapper<SamLine> toSamLines(String fname){
		ByteFile tf=ByteFile.makeByteFile(fname, false);
		LineParser1 lp=new LineParser1('\t');
		NameMapper<SamLine> list=new NameMapper<SamLine>();
		for(byte[] s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s[0]!='@'){
				SamLine sl=new SamLine(lp.set(s));
				list.put(new String(sl.rname()), sl);
			}
		}
		tf.close();
		return list;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Input reference sequence file. */
	private String in1=null;
	/** First primer SAM file. */
	private String sam1=null;
	/** Second primer SAM file. */
	private String sam2=null;
	/** Output sequence file. */
	private String out1=null;
	
	/** Emit one fake N read when no real interval is extracted. */
	private boolean ADD_FAKE_READS=true;
	/** Include primer spans rather than only the bases between them. */
	private boolean INCLUDE_PRIMERS=false;
	
	/** Input format descriptor. */
	private final FileFormat ffin1;
	/** Output format descriptor. */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	/** Maximum CRIS root records to request, or negative for unlimited. */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	/** Status and error stream redirected by the preparser when requested. */
	private java.io.PrintStream outstream=System.err;
	/** Emit verbose stream lifecycle messages. */
	public static boolean verbose=false;
	
}
