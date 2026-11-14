package synth;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import shared.Random;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Simulates multiple displacement amplification by repeatedly sampling random
 * fragments from a separator-delimited reference and appending them to an
 * amplified reference.  Optionally passes the amplified reference to
 * {@link RandomReads3} to synthesize error-containing reads.
 * <p>
 * Instances configure process-global FASTQ, compression and read-statistics
 * settings during construction, then run a single amplification pass.
 *
 * @author Brian Bushnell
 * @date Oct 17, 2014
 */
public class SynthMDA {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Runs one amplification pass and closes redirected status output.
	 * @param args input, reference-output and optional read-generation options
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		SynthMDA x=new SynthMDA(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Parses command-line options and creates input/output format descriptors.
	 * @param args command-line options; {@link PreParser} also expands config files
	 * @throws RuntimeException if no input reference is configured, or if the
	 *         amplified-reference output path cannot be written
	 */
	public SynthMDA(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		parser.build=7;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("minlen") || a.equals("ml")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("maxlen") || a.equals("mxl")){
				maxlen=Integer.parseInt(b);
			}else if(a.equals("cycles")){
				cycles=Integer.parseInt(b);
			}else if(a.equals("initialratio")){
				initialRatio=Float.parseFloat(b);
			}else if(a.equals("ratio")){
				ratio=Float.parseFloat(b);
			}else if(a.equals("refout")){
				out1=b;
			}else if(a.equals("perfect")){
				perfectrate=Float.parseFloat(b);
			}else if(a.equals("length")){
				readlength=Integer.parseInt(b);
			}else if(a.equals("paired")){
				paired=Parse.parseBoolean(b);
			}else if(a.equals("amp")){
				amp=Integer.parseInt(b);
			}
//			else if(a.equals("build")){
//				assert(false) : "Build should have been parsed by parser.";
//				build=Integer.parseInt(b);
//			}
			else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("prefix")){
				prefix=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
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
		minlen2=Tools.min(minlen2, minlen);
		
		{//Process parser fields
			Parser.processQuality();
			
			if(parser.maxReads>0){reads=parser.maxReads;}
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			if(ref==null){ref=parser.in1;}

			readsOut=parser.out1;
			
			extref=parser.extin;
			extout=parser.extout;
			build=parser.build;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(ref==null){throw new RuntimeException("Error - input reference must be specified.");}
		
		if(out1==null){
			out1=ReadWrite.stripToCore(ref)+"_"+Long.toHexString(shared.Shared.random().nextLong()&Long.MAX_VALUE)+".fa";
		}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffref=FileFormat.testInput(ref, FileFormat.FASTQ, extref, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Reads the input reference, amplifies it, writes the amplified reference and,
	 * when requested, invokes {@link RandomReads3} on that temporary reference.
	 * @param t timer to stop before printing throughput statistics
	 */
	void process(Timer t){
		
		ByteBuilder bb=new ByteBuilder();
		bb.append('$');
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ffref, null);
			if(verbose){outstream.println("Started cris");}
			cris.start();//4567
		}
		assert(!cris.paired());
		
		long readsProcessed=0;
		long basesProcessed=0;
		long maxContigLength=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffref==null || ffref.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					
					bb.append(r1.bases);
					bb.append('$');
					
					readsProcessed++;
					basesProcessed+=initialLength1;
					maxContigLength=Tools.max(maxContigLength, initialLength1);
				}
				
				final ArrayList<Read> listOut=reads;

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStream(cris);

		final long initialGoal=Tools.min((long)(bb.length()*initialRatio), 600000000L);
		//amplify() counts the separator it appends after each attempt, so
		//maxContigLength+1 is the longest possible accepted addition.
		if(initialGoal>0 && maxContigLength+1<minlen2){
			throw new IllegalArgumentException("At least one input contig must be long "+
					"enough to seed amplification: maxContigLength="+maxContigLength+
					", minlen2="+minlen2);
		}
		
		ByteBuilder dest=amplify(bb, false, minlen, maxlen, initialRatio);
		bb=null;
		for(int i=0; i<cycles; i++){
			dest=amplify(dest, i<1, minlen, maxlen, ratio);
//			if(dest.length()*ratio>1500000000){break;}
		}
//		assert(false) : cycles+", "+dest.length();
		
		TextStreamWriter tsw=(ffout1==null ? null : new TextStreamWriter(ffout1));
		if(tsw!=null){tsw.start();}
		
		bb=new ByteBuilder();
		for(int i=0, id=1; i<dest.length(); i++){
			byte b=dest.get(i);
			if(b=='$'){
				if(bb.length()>0){
					tsw.print(">"+id+"\n");
					tsw.println(bb.toString());
					id++;
				}
				bb.setLength(0);
			}else{
				bb.append(b);
			}
		}
		dest=null;
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(readsOut!=null){
			FileFormat ff=FileFormat.testOutput(readsOut, FileFormat.FASTQ, null, true, overwrite, false, false);
			assert(ff!=null);
			ArrayList<String> list=new ArrayList<String>();
			list.add("reads="+reads);
			list.add("length="+readlength);
			list.add("amp="+amp);
			if(paired){
				list.add("paired="+paired);
				list.add("interleaved="+paired);
			}
			list.add("build="+build);
			list.add("out="+readsOut);
			list.add("ow="+overwrite);
			list.add("minq="+16);
			list.add("midq="+25);
			list.add("maxq="+38);
			list.add("adderrors");
			list.add("snprate="+0.02);
			list.add("delrate="+0.005);
			list.add("insrate="+0.005);
			list.add("nrate="+0.005);
			list.add("maxinslen="+3);
			list.add("maxdellen="+3);
			list.add("maxnlen="+3);
			list.add("maxinss="+2);
			list.add("maxdels="+2);
			list.add("maxns="+2);
			list.add("maxsnps="+2);
			list.add("seed=-1");
			list.add("ref="+out1);
			if(prefix!=null){list.add("prefix="+prefix);}
			if(perfectrate>0){
				list.add("perfect="+perfectrate);
			}
			RandomReads3.main(list.toArray(new String[list.size()]));
		}
		
		boolean deleteRef=(readsOut!=null);
		if(deleteRef){
			if(verbose){System.err.println("Trying to delete "+out1);}
			try {
				File f=new File(out1);
				if(f.exists()){f.delete();}
			}catch(Exception e){
				// Best-effort cleanup; a leftover temporary reference is non-fatal.
				e.printStackTrace();
			}
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Extends a source buffer by sampling random forward or reverse-complement
	 * fragments until the requested ratio-derived goal is reached.
	 * @param source separator-delimited reference buffer to sample
	 * @param retain true to append to {@code source}; false to append to a fresh buffer
	 * @param minlen minimum attempted fragment length
	 * @param maxlen maximum attempted fragment length
	 * @param ratio multiplier used to derive the destination goal length
	 * @return the amplified buffer
	 */
	private ByteBuilder amplify(ByteBuilder source, boolean retain, int minlen, int maxlen, float ratio){
		assert(minlen<=maxlen && minlen>0 && maxlen>0);
		final int range=maxlen-minlen+1;
		final int slen=source.length();
		if(slen<minlen2*1.1f){
			KillSwitch.kill("Input ("+slen+") must be at least 10% longer than minlen ("+minlen2+").");
		}
		if(source.length()>=500000000){retain=false;}
		ByteBuilder dest=(retain ? source : new ByteBuilder());
		int goal=(int)Tools.min((long)(slen*ratio), 600000000);
		while(dest.length()<goal){
			final long initialLength=dest.length();
			final int start=randy.nextInt(slen);
			final int len0=minlen+randy.nextInt(range);
			final boolean forward=Tools.nextBoolean(randy);
			if(initialLength+(long)len0>1500000000){break;}
//			System.err.println(forward+", "+start+", "+len0);
			if(forward){
				final int stop=Tools.min(source.length(), start+len0);
//				System.err.println("stop="+stop);
				for(int i=start; i<stop; i++){
					byte b=source.get(i);
					if(b=='$'){
//						System.err.println("b="+(char)b);
						break;
					}
					dest.append(b);
				}
			}else{
				final int stop=Tools.max(0, start-len0);
//				System.err.println("stop="+stop);
				for(int i=start; i>=stop; i--){
					byte b=source.get(i);
					if(b=='$'){
//						System.err.println("b="+(char)b);
						break;
					}
					dest.append(AminoAcid.baseToComplementExtended[b]);
				}
			}
			dest.append('$');
			long added=dest.length()-initialLength;
//			System.err.println("added "+added+"/"+len0+" ("+initialLength+" -> "+dest.length()+")");
//			if(added<Tools.min(200, minlen) || (added<Tools.min(1000, minlen) && Tools.nextBoolean(randy))){dest.setLength(initialLength);}
			if(added<Tools.min(minlen2, minlen)){dest.setLength((int)initialLength);}
		}
		return dest;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String ref=null;
	private String out1=null;
	
	private String extref=null;
	private String extout=null;
	
	private final FileFormat ffref;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	private int minlen=10000;
	private int minlen2=4000;
	private int maxlen=150000;
	private int cycles=9;
	private float initialRatio=1.3f;
	private float ratio=2;
	
	private String prefix=null;
	
	/*--------------------------------------------------------------*/
	
	private long reads=12000000;
	private int readlength=150;
	private int amp=200;
	private boolean paired=true;
	private int build=7;
	private String readsOut=null;
	private float perfectrate=0;
	
	private final Random randy=Shared.threadLocalRandom();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
