package bloom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Aug 20, 2012
 *
 */
public class ErrorCorrect extends Thread{
	
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		String reads1=args[0];
		String reads2=(args.length>1 ? args[1] : null);
		
		int k=23;
		int cbits=4;
		int hashes=1;
		int thresh1=1;
		int thresh2=2;
		int matrixbits=34;
		long maxReads=-1;
		int buildpasses=1;
		long tablereads=-1; //How many reads to process when building the hashtable
		int buildStepsize=4;
		String output=null;
		boolean ordered=true;
		boolean overwrite=true;
		boolean append=false;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.equals("initialthresh") || a.equals("thresh1")){
				thresh1=Integer.parseInt(b);
			}else if(a.equals("thresh") || a.equals("thresh2")){
				thresh2=Integer.parseInt(b);
			}else if(a.startsWith("matrixbits")){
				matrixbits=Integer.parseInt(b);
			}else if(a.startsWith("hashes") || a.startsWith("multihash")){
				hashes=Integer.parseInt(b);
			}else if(a.startsWith("maxerrors")){
				ERROR_CORRECTION_LIMIT=Integer.parseInt(b);
			}else if(a.startsWith("passes")){
				buildpasses=Integer.parseInt(b);
			}else if(a.startsWith("stepsize") || a.startsWith("buildstepsize")){
				buildStepsize=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				System.err.println("Can't change threadcount for this class."); //THREADS=Integer.parseInt(b);
			}else if(a.startsWith("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.startsWith("tablereads")){
				tablereads=Parse.parseKMG(b);
			}else if(a.startsWith("build") || a.startsWith("genome")){
				Data.setGenome(Integer.parseInt(b));
				Data.sysout.println("Set genome to "+Data.GENOME_BUILD);
			}else if(a.equals("outputinfo") || a.startsWith("info")){
				OUTPUT_INFO=Parse.parseBoolean(b);
			}else if(a.startsWith("out")){
				output=b;
			}else if(a.startsWith("verbose")){
				KCountArray.verbose=Parse.parseBoolean(b);
//				verbose=KCountArray.verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		KCountArray kca=makeTable(reads1, reads2, k, cbits, hashes, buildpasses, matrixbits, tablereads, buildStepsize, thresh1, thresh2);
		
		detect(reads1, reads2, kca, k, thresh2, maxReads, output, ordered, append, overwrite);
	}
	
	//Only called by main...
	public static KCountArray makeTable(String reads1, String reads2, int k, int cbits, int hashes, int buildpasses, int matrixbits,
			long maxreads, int stepsize, int thresh1, int thresh2){
		
		Timer thash=new Timer();
		
		KmerCount6.maxReads=maxreads;
		int kbits=2*k;
		matrixbits=Tools.min(kbits, matrixbits);
		
		thash.start();
//		Data.sysout.println("kbits="+(kbits)+" -> "+(1L<<kbits)+", matrixbits="+(matrixbits)+" -> "+(1L<<matrixbits)+", cbits="+cbits+", gap="+gap+", hashes="+hashes);
		KCountArray kca=KCountArray.makeNew(1L<<matrixbits, cbits, hashes);
		
		if(buildpasses==1){
			
			KmerCount6.count(reads1, reads2, k, cbits, true, kca);
			kca.shutdown();
			
		}else{
			assert(buildpasses>1);
			KCountArray trusted=null;
			for(int i=1; i<buildpasses; i++){
				boolean conservative=i>2;// /*or, alternately, (trusted==null || trusted.capacity()>0.3)
				int step=(stepsize==1 ? 1 : stepsize+i%2);
//				if(!conservative){step=(step+3)/4;}
				if(!conservative){step=Tools.min(3, (step+3)/4);}
				
				KmerCount6.count(reads1, reads2, k, cbits, true, kca, trusted, maxreads, thresh1, step, conservative);
				
				kca.shutdown();
				Data.sysout.println("Trusted:   \t"+kca.toShortString());
				trusted=kca;
				kca=KCountArray.makeNew(1L<<matrixbits, cbits, hashes);
				
			}
			
			KmerCount6.count(reads1, reads2, k, cbits, true, kca, trusted, maxreads, thresh2, stepsize, true);
			
			kca.shutdown();
		}
		
		
		thash.stop();
//		Data.sysout.println(kca.toString());
		Data.sysout.println("Table :    \t"+kca.toShortString());
		Data.sysout.println("Hash time:      \t"+thash);
		return kca;
	}
	
	public static void detect(String reads1, String reads2, KCountArray kca, int k, int thresh, long maxReads, String output, boolean ordered, boolean append, boolean overwrite) {

		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
		}
		assert(cris!=null) : reads1;
		
		cris.start();
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		ConcurrentReadOutputStream ros=null;
		if(output!=null){
			String out1=output.replaceFirst("#", "1"), out2=null;
			
			if(cris.paired()){
				if(output.contains("#")){
					out2=output.replaceFirst("#", "2");
				}else{
					System.err.println("Writing interleaved.");
				}
			}
			
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, OUTPUT_INFO ? ".info" : null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, OUTPUT_INFO ? ".info" : null, true, overwrite, append, ordered);
			ros=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			
			assert(!ff1.sam()) : "Sam files need reference info for the header.";
		}
		
		
		if(ros!=null){
			ros.start();
			Data.sysout.println("Started output threads.");
		}
		
		detect(cris, kca, k, thresh, maxReads, ros);
		
		ReadWrite.closeStreams(cris, ros);
		if(verbose){System.err.println("Closed stream");}
	}
	
	public static void detect(ConcurrentReadInputStream cris, KCountArray kca, int k, int thresh, long maxReads, ConcurrentReadOutputStream ros) {
		Timer tdetect=new Timer();
		tdetect.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		long covered=0;
		long uncovered=0;
		
		long coveredFinal=0;
		long uncoveredFinal=0;
		
		long fullyCorrected=0;
		long failed=0;
		
		long totalBases=0;
		long totalReads=0;
		
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r : reads){
				Read r2=r.mate;
				{
					
//					if(r.numericID==23){verbose=true;}
					
					totalReads++;
					if(verbose){System.err.println();}
					totalBases+=r.length();
//					BitSet bs=detectErrors(r, kca, k, thresh);
					BitSet bs=detectErrorsBulk(r, kca, k, thresh, 1);
					if(verbose){System.err.println(toString(bs, r.length()));}
//					Data.sysout.println(toString(detectErrorsTips(r, kca, k, thresh), r.length()));
					if(verbose){System.err.println(toString(detectErrors(r, kca, k, thresh), r.length()-k+1));}
					if(bs==null){//No errors, or can't detect errors
						assert(false);
					}else{
						int x=bs.cardinality();
						covered+=x;
						uncovered+=(r.length()-x);
						if(x<r.length()){
							bs=correctErrors(r, kca, k, thresh, bs, ERROR_CORRECTION_LIMIT, MAX_ERROR_BURST);
						}
						int y=bs.cardinality();
						coveredFinal+=y;
						uncoveredFinal+=(r.length()-y);
						if(x<r.length()){
							if(y==r.length()){
								fullyCorrected++;
							}else{
								failed++;
							}
						}
					}
				}
				if(r2!=null){
					totalReads++;
					totalBases+=r2.length();
//					BitSet bs=detectErrors(r2, kca, k, thresh);
					BitSet bs=detectErrorsBulk(r2, kca, k, thresh, 1);
					if(verbose){System.err.println(toString(bs, r2.length()));}
//					Data.sysout.println(toString(detectErrorsTips(r2, kca, k, thresh), r2.length()));
					if(verbose){System.err.println(toString(detectErrors(r2, kca, k, thresh), r2.length()-k+1));}
					if(bs==null){//No errors, or can't detect errors
					}else{
						int x=bs.cardinality();
						covered+=x;
						uncovered+=(r2.length()-x);
						if(x<r2.length()){
							bs=correctErrors(r2, kca, k, thresh, bs, ERROR_CORRECTION_LIMIT, MAX_ERROR_BURST);
						}
						int y=bs.cardinality();
						coveredFinal+=y;
						uncoveredFinal+=(r2.length()-y);
						if(x<r2.length()){
							if(y==r2.length()){
								fullyCorrected++;
							}else{
								failed++;
							}
						}
					}
				}
			}
			
			if(ros!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
				if(DONT_OUTPUT_BAD_READS){removeBad(reads);}
				for(Read r : reads){
					if(r!=null){
						r.samline=null;
						assert(r.bases!=null);
						if(r.sites!=null && r.sites.isEmpty()){r.sites=null;}
					}
				}
//				System.err.println("Adding list of length "+readlist.size());
				ros.add(reads, ln.id);
			}
			
			cris.returnList(ln);
			//System.err.println("fetching list");
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(verbose){System.err.println("Finished reading");}
		cris.returnList(ln);
		if(verbose){System.err.println("Returned list");}
		
		tdetect.stop();
		Data.sysout.println("Detect time:    \t"+tdetect);
		Data.sysout.println("Total reads:    \t"+totalReads);
		Data.sysout.println("Total bases:    \t"+totalBases);
		Data.sysout.println("Reads Corrected:\t"+fullyCorrected);
		Data.sysout.println("Reads Failed:   \t"+failed);
		
		Data.sysout.println("\n - before correction - ");
		Data.sysout.println("Covered:        \t"+covered);
		Data.sysout.println("Uncovered:      \t"+uncovered);
		
		Data.sysout.println("\n -  after correction - ");
		Data.sysout.println("Covered:        \t"+coveredFinal);
		Data.sysout.println("Uncovered:      \t"+uncoveredFinal);
	}
	
	/** Sets a 1 bit at start of each kmer with count at least thresh */
	public static BitSet detectErrors(final Read r, final KCountArray kca, final int k, final int thresh){
		
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final byte[] bases=r.bases;
		
		int bslen=r.length()-k+1;
		if(bslen<1){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(bslen);
		
		int len=0;
		long kmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				
				if(len>=k){
					int count=kca.read(kmer);
					if(count>=thresh){
						bs.set(i+1-k);
					}
				}
			}
		}
		
		return bs;
	}
	
	/** Sets a 1 bit for every base covered by a kmer with count at least thresh */
	public static BitSet detectErrorsBulk(final Read r, final KCountArray kca, final int k, final int thresh, final int stepsize){
		
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final byte[] bases=r.bases;
		
		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.length());
		final int setlen=k;
		
		int len=0;
		long kmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				
				if(len>=k && ((len-k)%stepsize==0 || i==bases.length-1)){
					int count=kca.read(kmer);
					if(count>=thresh){
						bs.set(i+1-setlen, i+1);
					}
				}
			}
		}
		r.errors=bs.cardinality()-r.length();
		
		return bs;
	}
	
	/** Sets 1 for all bases.
	 * Then clears all bits covered by incorrect kmers. */
	public static BitSet detectTrusted(final Read r, final KCountArray kca, final int k, final int thresh, final int detectStepsize){
		
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final byte[] bases=r.bases;
		
		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.length());
		bs.set(0, r.length());
		final int setlen=k;
		
		int len=0;
		long kmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				
				if(len>=k && (i%detectStepsize==0 || i==bases.length-1)){
					int count=kca.read(kmer);
					if(count<thresh){
						bs.clear(i+1-setlen, i+1);
//						bs.clear(i+1-setlen+detectStepsize, i+1-detectStepsize);
//						bs.clear(i+k/2-detectStepsize, i+k/2+detectStepsize);
//						bs.clear(i+k/2);
					}
				}
			}
		}
//		assert(bases.length==r.length());
		return bs;
	}
	
	public static BitSet detectErrorsTips(final Read r, final KCountArray kca, final int k, final int thresh){
		
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final byte[] bases=r.bases;
		
		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.length());
		final int setlen=k;
		
		int len=0;
		long kmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				
				if(len>=k){
					int count=kca.read(kmer);
					if(count>=thresh){
						bs.set(i+1-setlen);
						bs.set(i);
					}
				}
			}
		}
		return bs;
	}
	

	/** Assumes bulk mode was used; e.g., any '0' bit is covered by no correct kmers */
	public static BitSet correctErrors(final Read r, final KCountArray kca, final int k, final int thresh, BitSet bs, final int maxCorrections, final int maxBurst){
		
		assert(!OUTPUT_INFO) : "TODO: Outputting correction data is not yet supported.";
		
		int corrections=0; //Alternately, corrections=r.errorsCorrected
		r.errors=0;
		
		if(bs.cardinality()==0){//Cannot be corrected
			r.errors=r.length();
			return bs;
		}

//		verbose=!bs.get(0);
		if(verbose){
			Data.sysout.println();
			Data.sysout.println(toString(bs, r.length()));
			Data.sysout.println(toString(detectErrorsTips(r, kca, k, thresh), r.length()));
			Data.sysout.println(toString(detectErrors(r, kca, k, thresh), r.length()-k+1));
		}
		
		
		int lastloc=-99;
		int burst=1;
		while(!bs.get(0) && corrections<maxCorrections){//While the read starts with a '0', correct from the right.
//			Data.sysout.println("Could not correct.");
//			return bs;
			int errorLoc=bs.nextSetBit(0)-1;//Location to left of first '1'
			if(Tools.absdif(errorLoc,lastloc)<=BURST_THRESH){burst++;}
			else{burst=1;}
			lastloc=errorLoc;
			boolean success=(burst<=MAX_ERROR_BURST) && correctFromRight(r, kca, k, thresh, bs, errorLoc);
			if(success){
				corrections++;
				bs=detectErrorsBulk(r, kca, k, thresh, 1);
				if(verbose){System.err.println(">\n"+toString(bs, r.length()));}
			}else{
				r.errors=r.length()-bs.cardinality();
//				r.errorsCorrected+=corrections;
				if(verbose){System.err.println("Could not correct.");}
				r.bases[errorLoc]='N';
				r.quality[errorLoc]=0;
				return bs;
			}
		}
		
		burst=1;
		while(bs.cardinality()<r.length() && corrections<maxCorrections){
			if(bs.get(0)){//First bit is a "1", can correct from the left
				int errorLoc=bs.nextClearBit(0);//Location to left of first '0'
				if(Tools.absdif(errorLoc,lastloc)<=BURST_THRESH){burst++;}
				else{burst=1;}
				lastloc=errorLoc;
				boolean success=(burst<=MAX_ERROR_BURST) && correctFromLeft(r, kca, k, thresh, bs, errorLoc);
				if(success){
					corrections++;
					bs=detectErrorsBulk(r, kca, k, thresh, 1);
					if(verbose){System.err.println(">\n"+toString(bs, r.length()));}
				}else{
					r.errors=r.length()-bs.cardinality();
//					r.errorsCorrected+=corrections;
					r.bases[errorLoc]='N';
					r.quality[errorLoc]=0;
					if(verbose){System.err.println("Could not correct.");}
					return bs;
				}
			}
		}

		r.errors=r.length()-bs.cardinality();
//		r.errorsCorrected+=corrections;
		assert(corrections<=maxCorrections);
		return bs;
	}
	
	
	/**
	 * @param r
	 * @param kca
	 * @param k
	 * @param thresh
	 * @param bs
	 * @param errorLoc
	 * @return
	 */
	private static boolean correctFromLeft(Read r, KCountArray kca, int k, int thresh, BitSet bs, int error) {
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final int setlen=k;
		final int startLoc=error-(setlen)+1;
		final byte oldBase=r.bases[error];
		final byte[] bases=r.bases;
		
		final int minAdvance=Tools.min(MIN_ADVANCE, bases.length-error);
		
		long kmer=0;
		int len=0;
		for(int i=startLoc; i<error; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				throw new RuntimeException("Can't correct from left!\nerror="+error+"\n"+toString(bs, bases.length)+"\n"+new String(bases)+"\nreadID: "+r.numericID);
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
			}
		}
		assert(len==setlen-1) : setlen+", "+len+", "+error+", "+startLoc;
		
		int[] counts=KillSwitch.allocInt1D(4);
		int[] dists=KillSwitch.allocInt1D(4);
		int maxLoc=Tools.min(bases.length-1, error+setlen-1);
		if(!bs.get(error+1)){maxLoc=Tools.min(maxLoc, error+9);}
		else{
			for(int i=error+2; i<=maxLoc; i++){
				if(!bs.get(i)){
					maxLoc=i-1;
					break;
				}
			}
		}
		
		if(verbose){System.err.println("correctFromLeft.  Error = "+error+", maxloc="+maxLoc);}
		for(int bnum=0; bnum<4; bnum++){
			byte c=AminoAcid.numberToBase[bnum];
			bases[error]=c;
			if(verbose){System.err.println("Considering "+(char)c);}
			long key=kmer;
			for(int loc=error; loc<=maxLoc; loc++){
				c=bases[loc];
				int x=AminoAcid.baseToNumber[c];
				if(x<0){
					if(verbose){System.err.println("break: N");}
					break;
				}
				key=((key<<2)|x)&mask;
				int count=kca.read(key);
				if(count<thresh){
					if(verbose){System.err.println("break: count="+count);}
					break;
				}
				dists[bnum]++;
				counts[bnum]+=count;
			}
		}
		bases[error]=oldBase;
		
		//Note:  I could require both to be the same, to decrease false-positives
		
		final int muid=maxUniqueIndex(dists);
		Arrays.sort(dists);
		final int advance=dists[3];
		final int delta=dists[3]-dists[2];
//		if(advance<minAdvance){return false;}
		if(delta<minAdvance){return false;}
		
		int best=(muid<0 ? maxUniqueIndex(counts) : muid);
		
		if(verbose){System.err.println("Best="+best+": "+Arrays.toString(dists)+"  \t"+Arrays.toString(counts));}
		if(best<0){return false;}
		byte bestC=AminoAcid.numberToBase[best];
		if(bestC==oldBase){return false;}
		bases[error]=bestC;
		
		r.quality[error]=(byte)Tools.min(10, 3+delta);
		
		return true;
	}
	
	
	
	/**
	 * @param r
	 * @param kca
	 * @param k
	 * @param thresh
	 * @param bs
	 * @param errorLoc
	 * @return
	 */
	private static boolean correctFromRight(Read r, KCountArray kca, int k, int thresh, BitSet bs, int error) {
		final int kbits=2*k;
		final int shift=kbits-2;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final int setlen=k;
		final int stopLoc=error+(setlen)-1;
		final byte oldBase=r.bases[error];
		final byte[] bases=r.bases;
		
		final int minAdvance=Tools.min(MIN_ADVANCE, error+1);
		
		long kmer=0;
		int len=0;
		
		for(int i=error+1; i<=stopLoc; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				throw new RuntimeException("Can't correct from right!\nerror="+error+"\n"+toString(bs, bases.length)+"\n"+new String(bases));
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
			}
//			Data.sysout.print((char)b);
		}
		kmer<<=2;
		
		if(verbose){
			Data.sysout.println();
			String s=Long.toBinaryString(kmer);
			while(s.length()<kbits){s="0"+s;}
			Data.sysout.println("kmer = \t"+s);
		}
		assert(len==setlen-1) : setlen+", "+len+", "+error+", "+stopLoc;
		
		int[] counts=KillSwitch.allocInt1D(4);
		int[] dists=KillSwitch.allocInt1D(4);
		int minLoc=Tools.max(0, error-setlen+1);
		if(error==0 || !bs.get(error-1)){minLoc=Tools.max(minLoc, error-9);}
		else{
			for(int i=error-2; i>=minLoc; i--){
				if(!bs.get(i)){
					minLoc=i+1;
					break;
				}
			}
		}
		
		if(verbose){
			Data.sysout.println("correctFromRight.  Error = "+error+", minloc="+minLoc);
			Data.sysout.println(new String(r.bases));
		}
		for(int bnum=0; bnum<4; bnum++){
			byte c=AminoAcid.numberToBase[bnum];
			bases[error]=c;
			if(verbose){System.err.println("Considering "+(char)c);}
			long key=kmer;
			for(int loc=error; loc>=minLoc; loc--){
				c=bases[loc];
				int x=AminoAcid.baseToNumber[c];
				if(x<0){
					if(verbose){System.err.println("break: N");}
					break;
				}
				key=((key>>2)|(((long)x)<<shift))&mask;
//				{
//					String s=Long.toBinaryString(key);
//					while(s.length()<kbits){s="0"+s;}
//					Data.sysout.println("mask="+Long.toBinaryString(mask)+", shift="+shift+", c="+c+", x="+x+", key  = \t"+s);
//				}
				int count=kca.read(key);
				if(count<thresh){
					if(verbose){System.err.println("break: count="+count);}
					break;
				}
				dists[bnum]++;
				counts[bnum]+=count;
			}
		}
		bases[error]=oldBase;
		
		//Note:  I could require both to be the same, to decrease false-positives
		
		final int muid=maxUniqueIndex(dists);
		Arrays.sort(dists);
		final int advance=dists[3];
		final int delta=dists[3]-dists[2];
//		if(advance<minAdvance){return false;}
		if(delta<minAdvance){return false;}
		
		int best=(muid<0 ? maxUniqueIndex(counts) : muid);
		
		if(verbose){System.err.println("Best="+best+": "+Arrays.toString(dists)+"  \t"+Arrays.toString(counts));}
		if(best<0){return false;}
		byte bestC=AminoAcid.numberToBase[best];
		if(bestC==oldBase){return false;}
		bases[error]=bestC;
		
		r.quality[error]=(byte)Tools.min(10, 3+delta);
		
		return true;
	}
	

	/** returns index of highest value, if unique; else a negative number */
	private static int maxUniqueIndex(int[] array){
		int max=array[0];
		int maxIndex=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>max){
				max=array[i];
				maxIndex=i;
			}else if(max==array[i]){
				maxIndex=-1;
			}
		}
		return maxIndex;
	}

	public static final String toString(BitSet bs, int len){
//		assert(verbose);
		StringBuilder sb=new StringBuilder(len);
		for(int i=0; i<len; i++){sb.append(bs.get(i) ? '1' : '0');}
		return sb.toString();
	}
	
	private static void removeBad(ArrayList<Read> list){
		
		for(int i=0; i<list.size(); i++){
			Read r=list.get(i);
			if(r.errors>0){
				if(r.mate==null || r.mate.errors>0){
					list.set(i, null);
				}
			}
		}
		
	}
	
	public static boolean verbose=false;
	/** Bails out if a read still has errors after correcting this many. */
	public static int ERROR_CORRECTION_LIMIT=6;
	/** Max allowed number of nearby corrections.
	 * A long error burst indicates the read simply has low coverage, and is not being corrected correctly. */
	public static int MAX_ERROR_BURST=3;
	/** Bursts have at most this distance between errors. E.G. '1' means errors are adjacent. */
	public static int BURST_THRESH=2;
	/** Withhold uncorrectable reads from output. */
	public static boolean DONT_OUTPUT_BAD_READS=false;
	/** Do not correct an error if it is at most this far from the next error.  Instead, bail out. */
	public static int MIN_ADVANCE=1;

	/** Number of threads used for error correction.  Does not control number of threads for creating the hash table.
	 * Additionally, up to 2 threads are used for reading and up to 2 for writing.  For this (singlethreaded) class, the number does nothing. */
	public static final int THREADS=1;
	
	/** Output correction data instead of the corrected read */
	public static boolean OUTPUT_INFO=false;
	
	
}
