package jgi;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import kmer.AbstractKmerTable;
import kmer.TableLoaderLockFree;
import kmer.TableReader;
import shared.Parse;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import template.BBTool_ST;

/**
 * @author Brian Bushnell
 * @date Mar 2, 2015
 *
 */
public class SplitNexteraLMP extends BBTool_ST {
	
	public static void main(String[] args){
		Timer t=new Timer();
		SplitNexteraLMP bbt=new SplitNexteraLMP(args);
		bbt.process(t);
	}
	
	@Override
	protected void setDefaults(){
		outStats="stderr";
		minReadLength=40;
		mask=false;
		merge=false;
		testmerge=0;
		pairedInput=true;
		symbol='J';
		useInnerLMP=false;
		RENAME=true;
	}

	/**
	 * @param args
	 */
	public SplitNexteraLMP(String[] args) {
		super(args);
		reparse(args);
		
		tables=(mask ? TableLoaderLockFree.makeTables(AbstractKmerTable.ARRAY1D, 12, -1L, false, 0.9) : null);
		
		if(outFrag1!=null && outFrag2==null && outFrag1.indexOf('#')>-1){
			outFrag2=outFrag1.replace("#", "2");
			outFrag1=outFrag1.replace("#", "1");
		}
		
		if(outUnk1!=null && outUnk2==null && outUnk1.indexOf('#')>-1){
			outUnk2=outUnk1.replace("#", "2");
			outUnk1=outUnk1.replace("#", "1");
		}
		
		if(testmerge>0){
			System.err.println("Testing merge rate.");
			float rate=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true).stdio() ? 1 : BBMerge.mergeableFraction(in1, in2, 1000000, 0.2f);
			merge=rate>0.1;
			System.err.println("Merge rate: "+Tools.format("%.2f%%", rate));
			if(!merge){
				System.err.println("Merging was disabled due to a low merge rate of "+Tools.format("%.3f", rate));
			}
		}
	}
	
	@Override
	public boolean parseArgument(String arg, String a, String b) {
		if(a.equals("symbol") || a.equals("junction")){
			assert(b!=null && b.length()==1) : "Junction symbol must be a single character.";
			symbol=(byte)b.charAt(0);
			return true;
		}else if(a.equals("outfrag") || a.equals("outfrag1") || a.equals("outf") || a.equals("outf1")){
			outFrag1=b;
			return true;
		}else if(a.equals("outfrag2") || a.equals("outf2")){
			outFrag2=b;
			return true;
		}else if(a.equals("outunknown") || a.equals("outunknown1") || a.equals("outu") || a.equals("outu1")){
			outUnk1=b;
			return true;
		}else if(a.equals("outunknown2") || a.equals("outu2")){
			outUnk2=b;
			return true;
		}else if(a.equals("outsingle") || a.equals("outs")){
			outSingle=b;
			return true;
		}else if(a.equals("minlen") || a.equals("minlength") || a.equals("ml")){
			minReadLength=Integer.parseInt(b);
			return true;
		}else if(a.equals("useinnerlmp") || a.equals("innerlmp")){
			useInnerLMP=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("rename")){
			RENAME=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("literal")){
			literals=(b==null ? null : b.split(","));
			return true;
		}else if(a.equals("mask")){
			mask=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("merge")){
			merge=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("testmerge")){
			testmerge=Double.parseDouble(b);
			if(testmerge>1){testmerge/=100;}
			return true;
		}else if(a.equals("rcomp")){
			rcomp=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("maskmiddle") || a.equals("mm")){
			maskMiddle=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("k")){
			k=Integer.parseInt(b);
			return true;
		}else if(a.equals("mink")){
			mink=Integer.parseInt(b);
			return true;
		}else if(a.equals("hdist") || a.equals("hammingdistance")){
			hdist=Integer.parseInt(b);
			return true;
		}else if(a.equals("hdist2") || a.equals("hammingdistance2")){
			hdist2=Integer.parseInt(b);
			return true;
		}else if(a.equals("edits") || a.equals("edist") || a.equals("editdistance")){
			edist=Integer.parseInt(b);
			return true;
		}else if(a.equals("dump")){
			dump=b;
			return true;
		}else if(a.equals("stats")){
			outStats=b;
			return true;
		}
		
		return false;
	}
	
	@Override
	protected void startupSubclass(){
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outFrag1, outFrag2, outUnk1, outUnk2)){
			throw new RuntimeException("\noverwrite="+overwrite+", append="+append+"\n" +
					"Can't write to output files "+out1+", "+out2+", "+outFrag1+", "+outFrag2+", "+outUnk1+", "+outUnk2+"\n");
		}
		
		if(!Tools.testForDuplicateFiles(true, in1, in2, qfin1, qfin2, out1, out2, qfout1, qfout2, outFrag1, outFrag2, outUnk1, outUnk2)){
			assert(false) : "Duplicate files.";
		}

		ffoutFrag1=FileFormat.testOutput(outFrag1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutFrag2=FileFormat.testOutput(outFrag2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffoutUnk1=FileFormat.testOutput(outUnk1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutUnk2=FileFormat.testOutput(outUnk2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffoutSingle=FileFormat.testOutput(outSingle, FileFormat.FASTQ, extout, true, overwrite, append, false);

		if(outFrag1!=null){
			final int buff=4;
			
			rosFrag=ConcurrentReadOutputStream.getStream(ffoutFrag1, ffoutFrag2, null, null, buff, null, false);
			rosFrag.start();
		}else{rosFrag=null;}
		
		if(outUnk1!=null){
			final int buff=4;
			
			rosUnk=ConcurrentReadOutputStream.getStream(ffoutUnk1, ffoutUnk2, null, null, buff, null, false);
			rosUnk.start();
		}else{rosUnk=null;}
		
		if(outSingle!=null){
			final int buff=4;
			
			rosSingle=ConcurrentReadOutputStream.getStream(ffoutSingle, null, null, null, buff, null, false);
			rosSingle.start();
		}else{rosSingle=null;}
		
	}
	
	
	@Override
	/** Iterate through the reads.
	 * This may optionally be overridden. */
	protected void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream rosLmp){
		
		if(mask){
			final TableLoaderLockFree loader=new TableLoaderLockFree(tables, k, mink, 0, hdist, edist, rcomp, maskMiddle);
			loader.setRefSkip(0);
			loader.hammingDistance2=hdist2;
			loader.storeMode(TableLoaderLockFree.SET_IF_NOT_PRESENT);
			long kmers=loader.processData(null, literals, false, false, false);
			outstream.println("Added "+kmers+" kmers.");
			if(dump!=null){
				ByteStreamWriter bsw=new ByteStreamWriter(dump, overwrite, false, true);
				bsw.start();
				for(AbstractKmerTable set : tables){
					set.dumpKmersAsBytes(bsw, k, 0, Integer.MAX_VALUE, null);
				}
				bsw.poisonAndWait();
			}
			
			reader=new TableReader(k, mink, 0, 0, 0, rcomp, maskMiddle);
			reader.trimSymbol=symbol;
			assert(kmers>0) : "There were no stored kmers; please check your settings.";
			
//			assert(false) : hdist+", "+hdist2+", "+maskMiddle+", "+rcomp+", "+k+", "+mink;
		}else{
			reader=null;
		}
		
		readsProcessed=0;
		basesProcessed=0;
		
		readsLmp=0;
		basesLmp=0;
		readsFrag=0;
		basesFrag=0;
		readsUnk=0;
		basesUnk=0;
		
		pairedInput=cris.paired();
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> listIn=(ln!=null ? ln.list : null);
			
			if(listIn!=null && !listIn.isEmpty()){
				Read r=listIn.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(listIn!=null && listIn.size()>0){
				if(verbose){outstream.println("Fetched "+listIn.size()+" reads.");}
				
				ArrayList<Read> outLmp=new ArrayList<Read>(listIn.size());
				ArrayList<Read> outFrag=new ArrayList<Read>(listIn.size());
				ArrayList<Read> outUnk=new ArrayList<Read>(listIn.size());
				ArrayList<Read> outSingle=new ArrayList<Read>(listIn.size());
				
				for(int idx=0; idx<listIn.size(); idx++){
					final Read r1=listIn.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					
					boolean keep=processReadPair(r1, r2, outLmp, outFrag, outUnk, outSingle);
					
				}

				for(Read r1 : outLmp){
					readsLmp+=r1.pairCount();
					basesLmp+=r1.pairLength();
					assert(r1.mate!=null);
					assert(r1.pairnum()==0);
					assert(r1.mate.pairnum()==1) : r1.mate.id+"\n\n"+r1.length()+"\n"+r1+"\n\n"+r1.mateLength()+"\n"+r1.mate+"\n\n";
					assert(r1.mate.mate==r1);
				}
				if(rosLmp!=null){rosLmp.add(outLmp, ln.id);}

				for(Read r1 : outFrag){
					readsFrag+=r1.pairCount();
					basesFrag+=r1.pairLength();
					assert(r1.mate!=null);
					assert(r1.pairnum()==0);
					assert(r1.mate.pairnum()==1);
					assert(r1.mate.mate==r1);
				}
				if(rosFrag!=null){rosFrag.add(outFrag, ln.id);}

				for(Read r1 : outUnk){
					readsUnk+=r1.pairCount();
					basesUnk+=r1.pairLength();
					assert(r1.mate!=null);
					assert(r1.pairnum()==0);
					assert(r1.mate.pairnum()==1);
					assert(r1.mate.mate==r1);
				}
				if(rosUnk!=null){rosUnk.add(outUnk, ln.id);}

				for(Read r1 : outSingle){
					readsSingle+=r1.pairCount();
					basesSingle+=r1.pairLength();
					assert(r1.pairnum()==0);
					assert(r1.mate==null);
				}
				if(rosSingle!=null){rosSingle.add(outSingle, ln.id);}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				listIn=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
	}
	
	@Override
	protected void showStatsSubclass(final Timer t, long readsIn, long basesIn){
		
		TextStreamWriter tsw=new TextStreamWriter(outStats==null ? "stderr" : outStats, overwrite, append, false);
		tsw.start();
		
		outstream.println("");
		
		final double rmult=(pairedInput ? 100.0 : 50.0)/readsIn;
		final double bmult=100.0/basesIn;
		
		//Note that this can go over 100%
		tsw.println("Long Mate Pairs:        \t"+readsLmp+" reads ("+Tools.format("%.2f",readsLmp*rmult)+"%) \t"+
				basesLmp+" bases ("+Tools.format("%.2f",basesLmp*bmult)+"%)");
		tsw.println("Fragment Pairs:         \t"+readsFrag+" reads ("+Tools.format("%.2f",readsFrag*rmult)+"%) \t"+
				basesFrag+" bases ("+Tools.format("%.2f",basesFrag*bmult)+"%)");
		tsw.println("Unknown Pairs:          \t"+readsUnk+" reads ("+Tools.format("%.2f",readsUnk*rmult)+"%) \t"+
				basesUnk+" bases ("+Tools.format("%.2f",basesUnk*bmult)+"%)");
		tsw.println("Singletons:             \t"+readsSingle+" reads ("+Tools.format("%.2f",readsSingle*100.0/readsIn)+"%) \t"+
				basesSingle+" bases ("+Tools.format("%.2f",basesSingle*bmult)+"%)");
		tsw.println("\n(Note: Read totals may exceed 100%, though base totals should not.)");
		tsw.println("");
		tsw.println("Adapters Detected:      \t"+junctionsDetected+" ("+Tools.format("%.2f%%)",junctionsDetected*100.0/junctionsSought));
		tsw.println("Bases Recovered:        \t"+(basesLmp+basesFrag+basesUnk+basesSingle)+
				" ("+Tools.format("%.2f%%)",(basesLmp+basesFrag+basesUnk+basesSingle)*bmult));
		if(merge){
			tsw.println("");
			tsw.println("Merged Pairs:           \t"+mergedReadCount+" ("+Tools.format("%.2f%%)",mergedReadCount*200.0/readsProcessed));
			tsw.println("Merged Bases:           \t"+mergedBaseCount+" ("+Tools.format("%.2f%%)",mergedBaseCount*100.0/basesProcessed));
		}
		
		errorState|=tsw.poisonAndWait();
	}
	
	@Override
	protected boolean processReadPair(Read r1, Read r2) {
		throw new RuntimeException("Do not use.");
	}
	
	boolean processReadPair(Read r1, Read r2, ArrayList<Read> outLmp, ArrayList<Read> outFrag, ArrayList<Read> outUnk, ArrayList<Read> outSingle) {
		boolean needsMasking=mask;
		if(merge){
			int insert=BBMerge.findOverlapStrict(r1, r2, false);
			if(insert>0){
				r2.reverseComplement();
				Read merged=r1.joinRead(insert);
				r2.reverseComplement();
				
				int a=1, b=0, c=0;
				if(mask){
					a=reader.kMask(merged, tables);
				}
				
				if(a>0 || true){
					mergedReadCount++;
					mergedBaseCount+=r1.length()+r2.length()-merged.length();
					return processMergedRead(merged, r1, r2, outLmp, outFrag, outUnk, outSingle);
				}else if(mask){
					needsMasking=false;
					b=reader.kMask(r1, tables);
					c=reader.kMask(r2, tables);
					if(b==0 && c==0){
						mergedReadCount++;
						mergedBaseCount+=r1.length()+r2.length()-merged.length();
						return processMergedRead(merged, r1, r2, outLmp, outFrag, outUnk, outSingle);
					}
				}
			}
		}
		
		if(needsMasking){
			int a=reader.kMask(r1, tables);
			int b=reader.kMask(r2, tables);
		}
		
		junctionsSought++;
		r1.start=Tools.indexOf(r1.bases, symbol);
		r1.stop=Tools.lastIndexOf(r1.bases, symbol);

		assert(r1==null || r1.pairnum()==0);
		assert(r2==null || r2.pairnum()==1);
		
		if(r2!=null){
			r2.start=Tools.indexOf(r2.bases, symbol);
			r2.stop=Tools.lastIndexOf(r2.bases, symbol);

			if(r1.start<0 && r2.start<0){
				if(verbose){System.err.println("Added unknown pair "+r1.id);}
				outUnk.add(r1);
				return true;
			}
			r1.mate=r2.mate=null;
		}else if(r1.start<0){
			if(verbose){System.err.println("Added singleton "+r1.id);}
			outSingle.add(r1);
			return true;
		}
		
		junctionsDetected++;
		
		Read r1left=null, r1right=null, r2left=null, r2right=null;
//		final Read r1left, r1right, r2left, r2right;
		
		if(r2==null){
			if(r1.start>=0){
				int left=r1.start;
				int right=r1.length()-r1.stop-1;
				
				r1left=(left>=minReadLength ? r1.subRead(0, r1.start) : null);
				r1right=null;
				r2left=null;
				r2right=(right>=minReadLength ? r1.subRead(r1.stop+1, r1.length()) : null);
				if(r2right!=null){
					r2right.setPairnum(1);
					if(RENAME){
						r2right.id=r2right.id.replaceFirst(" /1", " /2");
						r2right.id=r2right.id.replaceFirst(" 1:", " 2:");
					}
				}
			}
		}else if(r1.start>=0 && r2.start>=0){//confusing
			
			{
				int left=r1.start;
				int right=r1.length()-r1.stop-1;
				
				r1left=(left>=minReadLength ? r1.subRead(0, r1.start) : null);
				r1right=(right>=minReadLength ? r1.subRead(r1.stop+1, r1.length()) : null);
			}
			{
				int left=r2.start;
				int right=r2.length()-r2.stop-1;
				
				//Note these are reversed
				r2left=(right>=minReadLength ? r2.subRead(r2.stop+1, r2.length()) : null);
				r2right=(left>=minReadLength ? r2.subRead(0, r2.start) : null);
			}
		}else if(r1.start>=0){
			int left=r1.start;
			int right=r1.length()-r1.stop-1;
			
			r1left=(left>=minReadLength ? r1.subRead(0, r1.start) : null);
			r1right=(right>=minReadLength ? r1.subRead(r1.stop+1, r1.length()) : null);
			r2left=null;
			r2right=r2;
		}else if(r2.start>=0){
			int left=r2.start;
			int right=r2.length()-r2.stop-1;
			
			//Note these are reversed
			r2left=(right>=minReadLength ? r2.subRead(r2.stop+1, r2.length()) : null);
			r2right=(left>=minReadLength ? r2.subRead(0, r2.start) : null);
			r1left=r1;
			r1right=null;
		}else{
			assert(false) : r1.start+", "+r1.stop+(r2==null ? "null" : ", "+r2.start+", "+r2.stop);
		}
		
		boolean outerLMP=false, innerLMP=false, leftFrag=false, rightFrag=false;
		
		if(r1left!=null && r2right!=null){//outer lmp
			if(verbose){System.err.println("Added outer LMP "+r1.id);}
			r1left.mate=r2right;
			r2right.mate=r1left;
			outLmp.add(r1left);
			r1left=r2right=null;
			outerLMP=true;
		}
		
		if(r1right!=null && r2left!=null){//inner lmp
			if(verbose){System.err.println("Added inner LMP "+r1.id);}
			if(useInnerLMP){
				r1right.mate=r2left;
				r2left.mate=r1right;
				outLmp.add(r1right);
				r1right=r2left=null;
				innerLMP=true;
			}
		}
		
		if(r1left!=null && r2left!=null){//left frag
			if(verbose){System.err.println("Added left frag "+r1.id);}
			r1left.mate=r2left;
			r2left.mate=r1left;
			outFrag.add(r1left);
			r1left=r2left=null;
			leftFrag=true;
		}
		
		if(r1right!=null && r2right!=null){//right frag
			if(verbose){System.err.println("Added right frag "+r1.id);}
			r1right.mate=r2right;
			r2right.mate=r1right;
			outFrag.add(r1right);
			r1right=r2right=null;
			rightFrag=true;
		}
		
		//Singletons
		if(r1left!=null){
			if(verbose){System.err.println("Added singleton r1left "+r1left.id);}
			outSingle.add(r1left);
		}
		if(r1right!=null){
			if(verbose){System.err.println("Added singleton r1right "+r1right.id);}
			outSingle.add(r1right);
		}
		if(r2left!=null){
			if(verbose){System.err.println("Added singleton r2left "+r2left.id);}
			r2left.setPairnum(0);
			outSingle.add(r2left);
		}
		if(r2right!=null){
			if(verbose){System.err.println("Added singleton r2right "+r2right.id);}
			r2right.setPairnum(0);
			outSingle.add(r2right);
		}
		
		return true;
	}

	boolean processMergedRead(Read merged, Read r1, Read r2, ArrayList<Read> outLmp, ArrayList<Read> outFrag, ArrayList<Read> outUnk, ArrayList<Read> outSingle) {
		
//		int a=0, b, c;
//		if(mask){
//			a=reader.kMask(merged, tables);
//		}

		junctionsSought++;
		merged.start=Tools.indexOf(merged.bases, symbol);
		merged.stop=Tools.lastIndexOf(merged.bases, symbol);
		
		assert(merged!=null && merged.pairnum()==0);

		if(merged.start<0){
			if(verbose){System.err.println("Added frag "+r1.id);}
			outFrag.add(r1);
			return true;
		}

		junctionsDetected++;

		Read r1left=null, r1right=null;

		int left=merged.start;
		int right=merged.length()-merged.stop-1;

		r1left=(left>=minReadLength ? merged.subRead(0, merged.start) : null);
		r1right=(right>=minReadLength ? merged.subRead(merged.stop+1, merged.length()-1) : null);
		if(r1right!=null && r1left!=null){
			r1right.setPairnum(1);
			if(RENAME){
				r1right.id=r1right.id.replaceFirst(" /1", " /2");
				r1right.id=r1right.id.replaceFirst(" 1:", " 2:");
			}
			if(verbose){System.err.println("Added outer LMP "+merged.id);}
			r1left.mate=r1right;
			r1right.mate=r1left;
			outLmp.add(r1left);
			return true;
		}

		//Singletons
		if(r1left!=null){
			if(verbose){System.err.println("Added singleton r1left "+r1left.id);}
			outSingle.add(r1left);
		}
		if(r1right!=null){
			if(verbose){System.err.println("Added singleton r1right "+r1right.id);}
			outSingle.add(r1right);
		}

		return true;
	}
	
	@Override
	protected void shutdownSubclass(){
		errorState|=ReadWrite.closeStreams(null, rosFrag, rosUnk, rosSingle);
	}
	
	@Override
	protected final boolean useSharedHeader(){return true;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Masking Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	private String[] literals=new String[] {"CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG"};
	
	private final AbstractKmerTable[] tables;
	
	private int k=19;
	private int mink=11;
	private int hdist=1;
	private int hdist2=0;
	private int edist=0;
	private boolean rcomp=true;
	private boolean maskMiddle=false;
	
	private String dump=null;
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private TableReader reader;
	
	protected String outStats="stderr";
	
	protected String outFrag1;
	protected String outFrag2;
	
	protected String outUnk1;
	protected String outUnk2;
	
	protected String outSingle;
	
	protected FileFormat ffoutFrag1;
	protected FileFormat ffoutFrag2;
	
	protected FileFormat ffoutUnk1;
	protected FileFormat ffoutUnk2;
	
	protected FileFormat ffoutSingle;
	
	protected ConcurrentReadOutputStream rosFrag;
	protected ConcurrentReadOutputStream rosUnk;
	protected ConcurrentReadOutputStream rosSingle;

	private int minReadLength;

	private boolean mask;
	private boolean merge;
	private double testmerge;
	
	public long readsLmp(){return readsLmp;}
	public long basesLmp(){return basesLmp;}
	public long readsFrag(){return readsFrag;}
	public long basesFrag(){return basesFrag;}
	public long readsUnk(){return readsUnk;}
	public long basesUnk(){return basesUnk;}
	public long readsSingle(){return readsSingle;}
	public long basesSingle(){return basesSingle;}
	
	private long readsLmp=0;
	private long basesLmp=0;
	private long readsFrag=0;
	private long basesFrag=0;
	private long readsUnk=0;
	private long basesUnk=0;
	private long readsSingle=0;
	private long basesSingle=0;
	private long mergedReadCount=0;
	private long mergedBaseCount=0;
	
	private long junctionsSought=0, junctionsDetected=0;
	
	private boolean pairedInput;
	
	private byte symbol;
	private boolean useInnerLMP;
	
	private boolean RENAME;

}
