package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentHashMap;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray3;
import structures.IntList;
import structures.ListNum;
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * Masks a fasta file by inserting 'N' in place of low-complexity short repeats,
 * and anything covered by mapped reads in a sam file.
 * 
 * @author Brian Bushnell
 * @date Feb 18, 2014
 *
 */
public class BBMask{

	public static void main(String[] args){
		Timer t=new Timer();
		BBMask x=new BBMask(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	public BBMask(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(16);
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		boolean setEntropyMode=false;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.setThreads(b);
			}else if(a.equals("sampad") || a.equals("sampadding") || a.equals("sp")){
				samPad=Integer.parseInt(b);
			}else if(a.equals("entropymode")){
				entropyMode=Parse.parseBoolean(b);
				setEntropyMode=true;
			}else if(a.equals("maskrepeats") || a.equals("mr")){
				processRepeats=Parse.parseBoolean(b);
			}else if(a.equals("masklowentropy") || a.equals("masklowcomplexity") || a.equals("mlc") || a.equals("mle") || a.equals("me")){
				processEntropy=Parse.parseBoolean(b);
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1") || a.equals("ref")){
				inRef=b;
			}else if(a.equals("insam") || a.equals("samin") || a.equals("sam")){
				if(b!=null){
					for(String s : b.split(",")){
						inSam.add(s);
					}
				}
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1") || a.equals("output1")){
				outRef=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfinRef=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfoutRef=b;
			}else if(a.equals("extin")){
				extinRef=b;
			}else if(a.equals("extout")){
				extoutRef=b;
			}else if(a.equals("split")){
				splitMode=Parse.parseBoolean(b);
			}else if(a.equals("mink") || a.equals("kmin")){
				mink=mink2=Integer.parseInt(b);
			}else if(a.equals("maxk") || a.equals("kmax")){
				maxk=maxk2=Integer.parseInt(b);
			}else if(a.equals("k")){
				mink=maxk=mink2=maxk2=Integer.parseInt(b);
			}else if(a.equals("minkr") || a.equals("krmin")){
				mink=Integer.parseInt(b);
			}else if(a.equals("maxkr") || a.equals("krmax")){
				maxk=Integer.parseInt(b);
			}else if(a.equals("kr")){
				mink=maxk=Integer.parseInt(b);
			}else if(a.equals("mink2") || a.equals("kmin2") || a.equals("minke") || a.equals("kemin")){
				mink2=Integer.parseInt(b);
			}else if(a.equals("maxk2") || a.equals("kmax2") || a.equals("maxke") || a.equals("kemax")){
				maxk2=Integer.parseInt(b);
			}else if(a.equals("k2") || a.equals("ke")){
				mink2=maxk2=Integer.parseInt(b);
			}else if(a.equals("mincov")){
				mincov=Integer.parseInt(b);
			}else if(a.equals("maxcov")){
				maxcov=Integer.parseInt(b);
			}else if(a.equals("delcov") || a.equals("delcoverage")){
				includeDeletionCoverage=Parse.parseBoolean(b);
			}else if(a.equals("window") || a.equals("w")){
				window=Integer.parseInt(b);
			}else if(a.equals("ratio")){
				ratio=Float.parseFloat(b);
				if(!setEntropyMode){entropyMode=false;}
			}else if(a.equals("entropy") || a.equals("e")){
				entropyCutoff=Float.parseFloat(b);
				if(!setEntropyMode){entropyMode=true;}
			}else if(a.equals("lowercase") || a.equals("lc")){
				MaskByLowercase=Parse.parseBoolean(b);
			}else if(a.equals("minlen")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("mincount")){
				mincount=Integer.parseInt(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				stream.FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(!arg.contains("=") && FileFormat.hasSamOrBamExtension(arg)){
				inSam.add(arg);
			}
			
			else if(inRef==null && i==0 && Tools.looksLikeInputStream(arg)){
				inRef=arg;
			}
//			else if(outRef==null && i==1 && !arg.contains("=")){
//				outRef=arg;
//			}
			else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		bits32=(mincov>=Character.MAX_VALUE || maxcov>=Character.MAX_VALUE);
		
		{//Process parser fields
			Parser.processQuality();
		}

		assert(FastaReadInputStream.settingsOK());

		if(inRef==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(outRef!=null && outRef.equalsIgnoreCase("null")){outRef=null;}

		if(!Tools.testOutputFiles(overwrite, append, false, outRef)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+outRef+"\n");
		}
		
		ffoutRef=FileFormat.testOutput(outRef, FileFormat.FASTA, extoutRef, true, overwrite, append, false);

		ffinRef=FileFormat.testInput(inRef, FileFormat.FASTA, extinRef, true, true);
		
		if(inSam!=null && inSam.size()>0){
			ffinSam=new FileFormat[inSam.size()];
			for(int i=0; i<inSam.size(); i++){
				ffinSam[i]=FileFormat.testInput(inSam.get(i), FileFormat.SAM, null, true, false);
			}
		}else{
			ffinSam=null;
		}
		
		SamLine.CONVERT_CIGAR_TO_MATCH=false;
		
		EntropyTracker.verify=verify;
	}
	
	/*--------------------------------------------------------------*/


	public void process(Timer t0){
		
		Timer t=new Timer();
		{
			t.start();
			outstream.println("Loading input");
			
			boolean oldTI=FASTQ.TEST_INTERLEAVED;
			boolean oldFI=FASTQ.FORCE_INTERLEAVED;
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
			map=hashRef();
			FASTQ.TEST_INTERLEAVED=oldTI;
			FASTQ.FORCE_INTERLEAVED=oldFI;
			t.stop();
			
			outstream.println("Loading Time:                 \t"+t);
		}
		
		long repeats=0, lowcomplexity=0;
		long mapping=0;

		if(processRepeats && maxk>0){
			t.start();
			outstream.println("\nMasking repeats (to disable, set 'mr=f')");
//			repeats=maskRepeats_ST();
			repeats=maskRepeats();
			t.stop();
			
			double rpnano=refReads/(double)(t.elapsed);
			double bpnano=refBases/(double)(t.elapsed);
			
			String rpstring=""+refReads;
			String bpstring=""+refBases;
			String bmstring=""+repeats;
	
			while(rpstring.length()<12){rpstring=" "+rpstring;}
			while(bpstring.length()<12){bpstring=" "+bpstring;}
			while(bmstring.length()<12){bmstring=" "+bmstring;}
	
			outstream.println("Repeat Masking Time:          \t"+t);
			//outstream.println("Ref Scaffolds:          "+rpstring+" \t"+Tools.format("%.2fk scafs/sec", rpnano*1000000));
			outstream.println("Ref Bases:              "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
			outstream.println("Repeat Bases Masked:    "+bmstring);
		}

		if(processEntropy && maxk2>0){
			t.start();
			if(entropyMode){
				outstream.println("\nMasking low-entropy (to disable, set 'mle=f')");
				lowcomplexity=maskLowEntropy();
			}else{
				outstream.println("\nMasking low-complexity (to disable, set 'mlc=f')");
				lowcomplexity=maskLowComplexity(null);
			}
			t.stop();
			
			double rpnano=refReads/(double)(t.elapsed);
			double bpnano=refBases/(double)(t.elapsed);
			
			String rpstring=""+refReads;
			String bpstring=""+refBases;
			String bmstring=""+lowcomplexity;
	
			while(rpstring.length()<12){rpstring=" "+rpstring;}
			while(bpstring.length()<12){bpstring=" "+bpstring;}
			while(bmstring.length()<12){bmstring=" "+bmstring;}
	
			outstream.println("Low Complexity Masking Time:  \t"+t);
			//outstream.println("Ref Scaffolds:          "+rpstring+" \t"+Tools.format("%.2fk scafs/sec", rpnano*1000000));
			outstream.println("Ref Bases:              "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
			outstream.println("Low Complexity Bases:   "+bmstring);
		}
		
		if(ffinSam!=null){
			t.start();
			outstream.println("\nMasking from sam");
			mapping=maskSam();
			t.stop();
			
			double rpnano=samReads/(double)(t.elapsed);
			double bpnano=samBases/(double)(t.elapsed);
			
			String rpstring=""+samReads;
			String bpstring=""+samBases;
			String bmstring=""+mapping;
	
			while(rpstring.length()<12){rpstring=" "+rpstring;}
			while(bpstring.length()<12){bpstring=" "+bpstring;}
			while(bmstring.length()<12){bmstring=" "+bmstring;}
			
			outstream.println("Sam Masking Time:             \t"+t);
			outstream.println("Sam Reads Processed:    "+rpstring+" \t"+Tools.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Sam Bases Processed:    "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
			outstream.println("Sam Bases Masked:       "+bmstring);
		}
		long total=repeats+mapping+lowcomplexity, masked=0;
		
		if(total>0 || true){
			t.start();
			if(splitMode){
				masked=splitFromBitsets();
			}else{
				masked=maskFromBitsets(MaskByLowercase);
			}
			t.stop();
			outstream.println("Conversion Time:              \t"+t);
		}
		
		assert(total==masked) : repeats+", "+mapping+", "+lowcomplexity+", "+total+", "+masked;
		
		if(outRef!=null){
			t.start();
			outstream.println("\nWriting output");
			writeOutput();
			t.stop();
			outstream.println("Writing Time:                 \t"+t);
		}
		{
			t0.stop();
			String tstring=""+total;
			while(tstring.length()<12){tstring=" "+tstring;}
			outstream.println("\nTotal Bases Masked:     "+tstring+"/"+refBases+Tools.format("\t%.3f%%", total*100.0/refBases));
			outstream.println("Total Time:                   \t"+t0);
		}
		
		
		
		if(errorState){
			throw new RuntimeException("\nBBMask terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private static int setHighCoverage(BitSet bs, CoverageArray ca, int maxAllowedCoverage, int maxLen){
		int numSet=0;
		for(int i=0; i<maxLen; i++){
			if(ca.get(i)>maxAllowedCoverage){
				bs.set(i);
				numSet++;
			}
		}
		return numSet;
	}
	
	private static int setLowCoverage(BitSet bs, CoverageArray ca, int minAllowedCoverage, int maxLen){
		int numSet=0;
		for(int i=0; i<maxLen; i++){
			if(ca.get(i)<minAllowedCoverage){
				bs.set(i);
				numSet++;
			}
		}
		return numSet;
	}
	
	private long maskFromBitsets(final boolean lowercase){
		System.err.println("\nConverting masked bases to "+(lowercase ? "lower case" : "N")); //123
		long sum=0;
		for(Read r : map.values()){
			sum+=maskRead(r, (BitSet)r.obj, lowercase);
		}
		System.err.println("Done Masking");
		return sum;
	}
	
	public static int maskRead(final Read r, final BitSet bs, final boolean lowercase){
		final byte[] bases=r.bases;
		int sum=0;
		if(!lowercase){
			for(int i=0; i<bases.length; i++){
				if(bs.get(i)){
					if(bases[i]!='N'){sum++;}
					bases[i]='N';
				}else if(CONVERT_NON_ACGTN && !AminoAcid.isACGTN(bases[i])){
					bases[i]='N';
				}
			}
		}else{
			for(int i=0; i<bases.length; i++){
				if(bs.get(i)){
					if(!Tools.isLowerCase(bases[i]) && bases[i]!='N'){sum++;}
					bases[i]=(byte)Tools.toLowerCase(bases[i]);
				}else if(CONVERT_NON_ACGTN && !AminoAcid.isACGTN(bases[i])){
					bases[i]='N';
				}
			}
		}
		return sum;
	}
	
	private long splitFromBitsets(){
		System.err.println("\nSplitting reads by removing masked bases"); //123
		long sum=0;
		
		LinkedHashMap<String, Read> map2=new LinkedHashMap<String, Read>();

		for(String key : map.keySet()){
			Read r=map.get(key);
			BitSet bs=((BitSet)r.obj);
			int rnum=0;
			if(bs.isEmpty()){
				map2.put(key, r);
			}else{
				byte[] bases=r.bases;
				byte[] quals=r.quality;

				int lastGood, lastBad;
				if(bs.get(0)){
					lastGood=-1;
					lastBad=0;
				}else{
					lastGood=0;
					lastBad=-1;
				}
				
				int i=1;
				for(; i<bases.length; i++){
					if(bs.get(i)){
						if(lastGood==i-1){
							int len=lastGood-lastBad;
							if(len>0){
								byte[] bases2=KillSwitch.copyOfRange(bases, lastBad+1, i);
								byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, lastBad+1, i));
								Read r2=new Read(bases2, quals2, r.id+"_"+rnum, r.numericID, r.flags);
								Read old=map2.put(r2.id, r2);
								assert(old==null) : "Duplicate id "+r2.id; //TODO:  This can easily be resolved by making a new ID string.
							}
						}
						lastBad=i;
					}else{
						lastGood=i;
					}
				}
				if(lastGood==i-1){
					int len=lastGood-lastBad;
					if(len>0){
						byte[] bases2=KillSwitch.copyOfRange(bases, lastBad+1, i);
						byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, lastBad+1, i));
						Read r2=new Read(bases2, quals2, r.id+"_"+rnum, r.numericID, r.flags);
						Read old=map2.put(r2.id, r2);
						assert(old==null) : "Duplicate id "+r2.id; //TODO:  This can easily be resolved by making a new ID string.
					}
				}
			}
		}
		
		map.clear();
		map.putAll(map2);
		map2.clear();
		
		System.err.println("Done Splitting");
		return sum;
	}
	
	private void writeOutput(){
		
		ConcurrentReadOutputStream ros=null;
		if(ffoutRef!=null){
			final int buff=16;
			ros=ConcurrentReadOutputStream.getStream(ffoutRef, null, qfoutRef, null, buff, null, false);
			ros.start();
		}
		
		long i=0;
		for(String name : map.keySet()){
			Read r=map.get(name);
			ArrayList<Read> list=new ArrayList<Read>(1);
			list.add(r);
			if(ros!=null){ros.add(list, i);}
			i++;
		}
		errorState|=ReadWrite.closeStream(ros);
	}
	
	/*--------------------------------------------------------------*/
	
	private long maskSam(){
		long before=0, after=0;
		for(Read r : map.values()){
			before+=((BitSet)r.obj).cardinality();
		}
		for(FileFormat ff : ffinSam){
			//maskSam_ST(ff);
			maskSam_MT(ff);
		}
		for(Read r : map.values()){
			after+=((BitSet)r.obj).cardinality();
		}
		return after-before;
	}
	
	private void maskSam_MT(FileFormat ff){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null, null, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		
		MaskSamThread[] threads=new MaskSamThread[Shared.threads()];
//		outstream.println("Spawning "+numThreads+" threads.");
		for(int i=0; i<threads.length; i++){threads[i]=new MaskSamThread(cris);}
		for(int i=0; i<threads.length; i++){threads[i].start();}
		for(int i=0; i<threads.length; i++){
			while(threads[i].getState()!=Thread.State.TERMINATED){
				try {
					threads[i].join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		if(covmap!=null){//Move coverage array information to BitSets.
			for(String rs : map.keySet()){
				Read r=map.get(rs);
				BitSet bs=(BitSet)r.obj;
				CoverageArray ca=covmap.remove(rs);
				if(maxcov>-1){
					setHighCoverage(bs, ca, maxcov, r.length());
				}
				if(mincov>-1){
					setLowCoverage(bs, ca, mincov, r.length());
				}
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris);

		if(errorState){
			throw new RuntimeException("BBMask terminated in an error state; the output may be corrupt.");
		}
	}
	
	private class MaskSamThread extends Thread{
		
		MaskSamThread(ConcurrentReadInputStream cris_){
			cris=cris_;
		}
		
		@Override
		public void run(){
			maskSam(cris);
		}
		
		final ConcurrentReadInputStream cris;
		
		
	}
	
	private void maskSam_ST(FileFormat ff){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null, null, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		
		maskSam(cris);

		errorState|=ReadWrite.closeStreams(cris);

		if(errorState){
			throw new RuntimeException("BBMask terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void maskSam(ConcurrentReadInputStream cris){
		
		long samReads=0;
		long samBases=0;
		IntList ranges=new IntList(16);
		
		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					assert(r.mate==null);

					final int initialLength1=r.length();

					{
						samReads++;
						samBases+=initialLength1;
						
						if(r.mapped()){
							SamLine sl=r.samline;
							assert(sl!=null) : "No sam line for read "+r;
							byte[] rname=sl.rname();
							assert(rname!=null) : "No rname for sam line "+sl;
							final String rs=new String(rname);
							Read ref=map.get(rs);
							if(ref==null){
								handleNoRef(rs);
							}else{
								final int reflen=ref.length();
								assert(ref!=null) : "Could not find reference scaffold '"+rs+"' for samline \n"+sl+"\n in set \n"+map.keySet()+"\n";
								if(covmap==null){
									BitSet bs=(BitSet)ref.obj;
									mask(bs, sl, reflen);
								}else{
									CoverageArray ca=covmap.get(rs);
									increment(ca, sl, r.match, reflen, ranges, includeDeletionCoverage, samPad);
								}
							}
						}
					}
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		synchronized(this){
			this.samBases+=samBases;
			this.samReads+=samReads;
		}
	}
	
	private void mask(BitSet bs, SamLine sl, int reflen){
		final int start=Tools.max(0, sl.start(true, false)-samPad);
		final int stop=Tools.min(sl.stop(start, true, false)+1+samPad, reflen);
		if(stop>start){
			synchronized(bs){//Potential bottleneck
				bs.set(start, stop);
			}
		}
	}
	
	public static void increment(CoverageArray ca, SamLine sl, byte[] match, int reflen, IntList ranges, boolean includeDels, int samPad){
		final int start=Tools.max(0, sl.start(true, false));
		final int stop=Tools.min(sl.stop(start, true, false)+1, reflen);
		if(stop>start){
			ranges.clear();
			boolean numeric=false;
			if(match==null){
				assert(sl.cigar!=null);
				match=sl.toShortMatch(true);
				numeric=true;
			}else{
				for(byte b : match){
					if(Tools.isDigit(b)){
						numeric=true;
						break;
					}
				}
			}
			if(numeric){match=Read.toLongMatchString(match);}
			fillRanges(match, start, stop, ranges, includeDels);
//			assert(false) : ranges;
			if(ranges.size>0){
				if(samPad!=0){//Pad the ranges, but don't let them overlap
					ranges.set(0, Tools.mid(0, reflen, ranges.get(0)-samPad));
					ranges.set(1, Tools.mid(0, reflen, ranges.get(1)+samPad));
					for(int i=2; i<ranges.size; i+=2){
						ranges.set(i, Tools.mid(ranges.get(i-1), reflen, ranges.get(i)-samPad));
						ranges.set(i+1, Tools.mid(0, reflen, ranges.get(i+1)+samPad));
					}
				}
				synchronized(ca){//Potential bottleneck
					ca.incrementRanges(ranges, 1);
				}
			}
		}
	}
	
	public static void fillRanges(byte[] longmatch, int start, int stop, IntList ranges, boolean includeDels){
		assert(ranges.size==0);
		byte mode='?', lastMode='?';
		int rpos=start;
		int lastRpos=start;
		int rstart=start;
		for(int mpos=0; mpos<longmatch.length; mpos++){
			byte m=longmatch[mpos];
			if(m=='m' || m=='s' || m=='S' || m=='N' || m=='B'){//Little 's' is for a match classified as a sub to improve the affine score.
				mode='m';
				rpos++;
			}else if(m=='I' || m=='X' || m=='Y'){
				mode='I';
			}else if(m=='D'){
				mode='D';
				rpos++;
			}else if(m=='C'){
				mode='C';
				rpos++;
			}else{
				throw new RuntimeException("Invalid match string character '"+(char)m+"' = "+m+" (ascii).  " +
						"Match string should be in long format here.");
			}
			if(mode!=lastMode){
				if(mpos>0){
					if(lastMode=='m'){
						ranges.add(rstart);
						ranges.add(lastRpos);
					}else if(mode=='D' && includeDels){
						ranges.add(rstart);
						ranges.add(lastRpos);
					}
				}
				rstart=lastRpos;
			}
			lastMode=mode;
			lastRpos=rpos;
		}
//		assert(false) : rstart+", "+rpos+", "+(char)mode+", "+(char)lastMode;
		//Final cycle
		if(lastMode=='m'){
			ranges.add(rstart);
			ranges.add(lastRpos);
		}else if(mode=='D' && includeDels){
			ranges.add(rstart);
			ranges.add(lastRpos);
		}
		assert(rpos==stop) : start+", "+stop+", "+rpos+", "+rstart+", "+lastRpos+"\n"+new String(longmatch)+"\n"+ranges;
		assert((ranges.size&1)==0);
	}
	
	private void handleNoRef(String rname){
		assert(rname!=null);
		String ret=norefSet.putIfAbsent(rname, rname);
		if(ret==null){
			System.err.println("Warning! Scaffold not found in assembly: "+rname);
		}
	}
	
	/*--------------------------------------------------------------*/

	private LinkedHashMap<String, Read> hashRef(){


		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffinRef, null, qfinRef, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		
		final LinkedHashMap<String, Read> hmr=new LinkedHashMap<String, Read>();
		if(ffinSam!=null && (mincov>=0 || maxcov>=0)){
			covmap=new HashMap<String, CoverageArray>();
		}
		
		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					final int len=r.length();
					final byte[] bases=r.bases;
					
//					refLengths.add(len);
					final BitSet bs=new BitSet(len);
					r.obj=bs;
					
					if(covmap!=null){
						if(bits32){
							covmap.put(r.id, new CoverageArray3(covmap.size(), len));
						}else{
							covmap.put(r.id, new CoverageArray2(covmap.size(), len));
						}
					}
					
					if(MaskByLowercase){
						for(int i=0; i<len; i++){
							if(bases[i]=='N' || Tools.isLowerCase(bases[i])){bs.set(i);}
						}
					}else{
						for(int i=0; i<len; i++){
							if(bases[i]=='N'){bs.set(i);}
						}
					}
					
					refReads++;
					refBases+=len;
					Read old=hmr.put(r.id, r);
					assert(old==null) : "Duplicate reference scaffold name "+r.id;
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		errorState|=ReadWrite.closeStreams(cris);

		if(errorState){
			throw new RuntimeException("BBMask terminated in an error state; the output may be corrupt.");
		}
		
		return hmr;
	}
	
	/*--------------------------------------------------------------*/
	
	
	private long maskLowComplexity(short[][] matrix){
		long sum=0;
		if(matrix==null){matrix=new short[16][];}
		for(Read r : map.values()){
			sum+=maskLowComplexity(r, mink2, maxk2, window, ratio, matrix);
		}
		return sum;
	}
	
	private static int maskLowComplexity(Read r, int mink, int maxk, int window, float ratio, short[][] matrix){
		
		final byte[] bases=r.bases;
		final BitSet bs=(BitSet)r.obj;
		
		int before=bs.cardinality();
//		System.err.println("\nbefore="+before+"\n"+new String(bases)+"\n"+bs);
		
		for(int k=mink; k<=maxk; k++){
			if(matrix[k]==null){matrix[k]=new short[(1<<(2*k))];}
		}
		
		for(int k=mink; k<=maxk; k++){
			final short[] counts=matrix[k];
			final int kmerspace=(1<<(2*k));
			final int mincount=(int)Math.ceil(ratio*Tools.min(window, kmerspace));
			maskLowComplexity(bases, bs, k, window, mincount, counts);
		}

		int after=bs.cardinality();
		
//		System.err.println("before="+before+", after="+after+"\n"+new String(bases)+"\n"+bs);
		
		return after-before;
	}
	
	
	private static void maskLowComplexity(final byte[] bases, final BitSet bs, final int k, final int window, final int mincount, final short[] counts){
		assert(k>0) : "k must be greater than 0";
		
		if(verify){
			for(int c : counts){assert(c==0);}
		}
		
		final int mask=(k>15 ? -1 : ~((-1)<<(2*k)));
		int current=0, ns=0;
		int kmer=0, kmer2=0;
		
		for(int i=0, i2=-window; i2<bases.length; i++, i2++){
			
//			System.err.println("\nStart: i="+i+", current="+current+", ns="+ns+"\n"+Arrays.toString(counts));
			
			if(i<bases.length){
				final byte b=bases[i];
				final int n=Dedupe.baseToNumber[b];
				kmer=((kmer<<2)|n)&mask;
				
				if(!AminoAcid.isFullyDefined(b)){ns++;}
				if(counts[kmer]<1){
					assert(counts[kmer]==0);
					current++;
				}
				counts[kmer]++;
				if(verify){assert(current==Tools.cardinality(counts)) : current+", "+Tools.cardinality(counts)+"\n"+Arrays.toString(counts);}
				
//				System.err.println("Added "+kmer+"; counts["+kmer+"]="+counts[kmer]);
			}
			
			if(i2>=0){
				final byte b2=bases[i2];
				final int n2=Dedupe.baseToNumber[b2];
				kmer2=((kmer2<<2)|n2)&mask;

				if(!AminoAcid.isFullyDefined(b2)){
					ns--;
					assert(ns>=0);
				}
				counts[kmer2]--;
				if(counts[kmer2]<1){
					assert(counts[kmer2]==0) : Arrays.toString(counts);
					current--;
				}
				if(verify){assert(current==Tools.cardinality(counts)) : current+", "+Tools.cardinality(counts)+"\n"+Arrays.toString(counts);}
				
//				System.err.println("Removed "+kmer2+"; count="+counts[kmer2]);
			}
			
			if(verify && i2>-1 && i<bases.length){
				assert(shared.Vector.sum(counts)==window);
			}
			
			if(current<mincount && ns<1 && i2>=-1 && i<bases.length){
//				System.err.println("Masked ("+(i2+1)+", "+(i+1)+")");
				bs.set(i2+1, i+1);
			}
		}
		
	}
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	
	private long maskLowEntropy(){
		ArrayBlockingQueue<Read> queue=new ArrayBlockingQueue<Read>(map.size());
		for(Read r : map.values()){queue.add(r);}
		int numThreads=Tools.min(Shared.threads(), queue.size());
		MaskLowEntropyThread[] threads=new MaskLowEntropyThread[numThreads];
		long sum=0;
//		outstream.println("Spawning "+numThreads+" threads.");
		for(int i=0; i<threads.length; i++){threads[i]=new MaskLowEntropyThread(queue, mink2, maxk2, window, entropyCutoff);}
		for(int i=0; i<threads.length; i++){threads[i].start();}
		for(int i=0; i<threads.length; i++){
			while(threads[i].getState()!=Thread.State.TERMINATED){
				try {
					threads[i].join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			sum+=threads[i].masked;
		}
		return sum;
	}
	
	private class MaskLowEntropyThread extends Thread{
		
		MaskLowEntropyThread(ArrayBlockingQueue<Read> queue_, int mink_, int maxk_, int window_, float cutoff_){
			queue=queue_;
			minkT=mink_;
			maxkT=maxk_;
			windowT=(short)window_;
			cutoff=cutoff_;
			trackers=new EntropyTracker[maxkT+1];
			for(int i=minkT; i<=maxkT; i++){
				trackers[i]=new EntropyTracker(i, windowT, false, cutoff, true);
			}
		}
		
		@Override
		public void run(){
			for(Read r=queue.poll(); r!=null; r=queue.poll()){
				masked+=maskLowEntropy(r, minkT, maxkT, trackers);
			}
		}
		
		final ArrayBlockingQueue<Read> queue;
		final int minkT;
		final int maxkT;
		final float cutoff;
		final int windowT;
		final EntropyTracker[] trackers;
		long masked=0;
		
	}
	
	private int maskLowEntropy(Read r, int mink, int maxk, EntropyTracker[] trackers){
//		outstream.println("maskLowEntropy("+r.numericID+", "+mink+", "+maxk+", "+window+", "+cutoff+", "+matrix.length+", "+countCounts.length+")");
//		System.err.println(new String(r.bases));
		final byte[] bases=r.bases;
		final BitSet bs=(BitSet)r.obj;
		
		int before=bs.cardinality();
		
		for(int k=mink; k<=maxk; k++){
			maskLowEntropy(bases, bs, trackers[k]);
		}

		int after=bs.cardinality();
		
		return after-before;
	}
	
	
	private static void maskLowEntropy(final byte[] bases, final BitSet bs, final EntropyTracker et){
		final int window=et.windowBases();
		if(bases.length<window){return;}
		
		et.clear();
		for(int i=0, min=window-1; i<bases.length; i++){
			et.add(bases[i]);
			if(i>=min && et.ns()<1 && !et.passes()){bs.set(et.leftPos(), et.rightPos()+1);}
		}
	}
	
	
	/*--------------------------------------------------------------*/
	


	private long maskRepeats_ST(){
		long sum=0;
		for(Read r : map.values()){
			sum+=maskRepeats(r, mink, maxk, mincount, minlen);
		}
		return sum;
	}
	
	private long maskRepeats(){
		ArrayBlockingQueue<Read> queue=new ArrayBlockingQueue<Read>(map.size());
		for(Read r : map.values()){queue.add(r);}
		int numThreads=Tools.min(Shared.threads(), queue.size());
		MaskRepeatThread[] threads=new MaskRepeatThread[numThreads];
		long sum=0;
		for(int i=0; i<threads.length; i++){threads[i]=new MaskRepeatThread(queue, mink, maxk, mincount, minlen);}
		for(int i=0; i<threads.length; i++){threads[i].start();}
		for(int i=0; i<threads.length; i++){
			while(threads[i].getState()!=Thread.State.TERMINATED){
				try {
					threads[i].join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			sum+=threads[i].masked;
		}
		return sum;
	}
	
	private class MaskRepeatThread extends Thread{
		
		MaskRepeatThread(ArrayBlockingQueue<Read> queue_, int mink_, int maxk_, int mincount_, int minlen_){
			queue=queue_;
			minkT=mink_;
			maxkT=maxk_;
			mincountT=mincount_;
			minlenT=minlen_;
		}
		
		@Override
		public void run(){
			for(Read r=queue.poll(); r!=null; r=queue.poll()){
				masked+=maskRepeats(r, minkT, maxkT, mincountT, minlenT);
			}
		}
		
		final ArrayBlockingQueue<Read> queue;
		final int minkT;
		final int maxkT;
		final int mincountT;
		final int minlenT;
		long masked=0;
		
	}
	
	private static int maskRepeats(Read r, int mink, int maxk, int mincount, int minlen){
		final byte[] bases=r.bases;
		final BitSet bs=(BitSet)r.obj;
		
		int before=bs.cardinality();
//		System.err.println("\nbefore="+before+"\n"+new String(bases)+"\n"+bs);
		
		for(int k=mink; k<=maxk; k++){
			maskRepeats(bases, bs, k, Tools.max(minlen, k*mincount));
		}

		int after=bs.cardinality();
		
//		System.err.println("before="+before+", after="+after+"\n"+new String(bases)+"\n"+bs);
		
		return after-before;
	}
	
	
	private static void maskRepeats(final byte[] bases, final BitSet bs, final int k, final int minlen){
		final int lim=bases.length-k;
		final int mask=(k>15 ? -1 : ~((-1)<<(2*k)));
		for(int loc=0; loc<lim; loc++){
			int len=repeatLength(bases, k, mask, loc);
			if(len>=minlen){
				int a=loc-k, b=loc-k+len;
				bs.set(a, b);
//				System.err.println("len="+len+", minlen="+minlen+", set "+(loc-k)+"-"+(loc-k+len));
				loc=Tools.max(loc, b-minlen);
//				System.err.println("Reset loc to "+loc);
			}else{
//				System.err.println("len="+len+" < minlen="+minlen);
			}
		}
		
	}
	
	
	private static int repeatLength(final byte[] bases, final int k, final int mask, final int loc){
		
		final int lim=bases.length;
		final int key=getInitialKey(bases, loc, k);
		if(key<0){return 0;}
		int kmer=key;
		int gap=0, last=-1;
		for(int i=loc; i<lim && gap<k; i++){
			final byte b=bases[i];
			final int n=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|n)&mask;
			if(kmer==key){
				last=i;
				gap=0;
			}else{
				gap++;
			}
//			System.err.println("i="+i+", lim="+lim+", gap="+gap+", last="+last+", b="+(char)b+", n="+n+", key="+key+", kmer="+kmer);
		}
		
//		System.err.println("k="+k+", mask="+mask+", loc="+loc+", last="+last);
		
		return (last<0 ? 0 : last-loc+k+1);
	}
	
	private static int getInitialKey(byte[] bases, int loc, int k){
		assert(k<16);
		int start=loc-k;
		int key=0;
		if(start<0){return -1;}
		for(int i=start; i<loc; i++){
			final byte b=bases[i];
			final int n=Dedupe.baseToNumber[b];
			key=(key<<2)|n;
		}
		assert(key>=0);
		return key;
	}

	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/

	private LinkedHashMap<String, Read> map=null;
	private ConcurrentHashMap<String, String> norefSet=new ConcurrentHashMap<String, String>(256, .75f, 16);
	private HashMap<String, CoverageArray> covmap=null;
//	private IntList refLengths=new IntList();
	
	private long refReads=0;
	private long refBases=0;
//	private long repeatsMasked=0;
	
	private long samReads=0;
	private long samBases=0;
//	private long samMasked=0;
	
	public boolean errorState=false;

	private String inRef=null;
	private ArrayList<String> inSam=new ArrayList<String>();

	private String qfinRef=null;

	private String outRef=null;

	private String qfoutRef=null;

	private String extinRef=null;
	private String extoutRef=null;
	
	private boolean overwrite=true;
	private boolean append=false;

	private long maxReads=-1;

	private boolean processRepeats=false;
	private int mink=5;
	private int maxk=5;
	private int minlen=40;
	private int mincount=4;

	private boolean processEntropy=true;
	private boolean entropyMode=true;
	private boolean splitMode=false;
	private int mink2=5;
	private int maxk2=5;
	private int window=80;
	private float ratio=0.35f; //For complexity, if not in entropyMode
	private float entropyCutoff=0.70f;
	
	/** Use 32-bit coverage arrays */
	private boolean bits32=true;
	/** Include deletions when calculating coverage */
	private boolean includeDeletionCoverage=true;

	/** If nonnegative, mask bases with coverage outside this range. */
	private int mincov=-1;
	/** If nonnegative, mask bases with coverage outside this range. */
	private int maxcov=-1;
	
	private int samPad=0;

	private final FileFormat ffinRef;
	private final FileFormat[] ffinSam;

	private final FileFormat ffoutRef;

	private PrintStream outstream=System.err;

	/*--------------------------------------------------------------*/
	
	public static boolean verbose=false;
	public static boolean CONVERT_NON_ACGTN=true;
	private static boolean verify=false;
	private static boolean MaskByLowercase=false;

}
