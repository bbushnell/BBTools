package var;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import align2.MultiStateAligner9ts;
import align2.TranslateColorspaceRead;
import dna.Data;
import fileIO.ReadWrite;
import fileIO.TextFile;
import pacbio.CalcCoverageFromSites;
import pacbio.SiteR;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentLegacyReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;
import stream.SiteScoreR;
import structures.ListNum;

/** Splits output files across blocks for low memory usage */
public class GenerateVarlets2 {
	
	
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Data.GENOME_BUILD=-1;
		
		String reads1=args[0];
		String reads2=args[1].equalsIgnoreCase("null") ?  null : args[1];
		String outname=args[2];
//		assert(outname.contains("#"));
		
		String sitesfile=null;
		
		byte minChrom=1;
		byte maxChrom=1;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		
		
		for(int i=3; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(a.equals("condense")){
				CONDENSE=Parse.parseBoolean(b);
			}else if(a.equals("condensesnps")){
				CONDENSE_SNPS=Parse.parseBoolean(b);
			}else if(a.startsWith("splitsubs")){
				SPLIT_SUBS=Parse.parseBoolean(b);
			}else if(a.equals("tosssolo1")){
				TOSS_SOLO1=Parse.parseBoolean(b);
			}else if(a.equals("tosssolo2")){
				TOSS_SOLO2=Parse.parseBoolean(b);
			}else if(a.startsWith("minchrom")){
				minChrom=Byte.parseByte(b);
			}else if(a.startsWith("maxchrom")){
				maxChrom=Byte.parseByte(b);
			}else if(a.startsWith("build") || a.startsWith("genomebuild") || a.startsWith("genome")){
				Data.setGenome(Integer.parseInt(b));
				System.out.println("Set GENOME_BUILD to "+Data.GENOME_BUILD);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(Integer.parseInt(b));
			}else if(a.startsWith("buffer") || a.startsWith("writebuffer")){
				WRITE_BUFFER=(Integer.parseInt(b));
			}else if(a.startsWith("maxreads")){
				MAX_READS=(Long.parseLong(b));
			}else if(a.equals("blocksize")){
				BLOCKSIZE=(Integer.parseInt(b));
			}else if(a.startsWith("sites") || a.startsWith("sitesfile")){
				sitesfile=b;
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		assert(minChrom<=maxChrom && minChrom>=0);
		if(Data.GENOME_BUILD<0){throw new RuntimeException("Please set genome number.");}
		
		if(ReadWrite.ZIPLEVEL<2){ReadWrite.ZIPLEVEL=2;}
		GenerateVarlets2 gv=new GenerateVarlets2(reads1, reads2, outname, minChrom, maxChrom, MAX_READS, sitesfile);
		gv.process();
	}
	
	public GenerateVarlets2(String fname1, String fname2, String outname_, byte minChrom, byte maxChrom, long maxReads, String sitesfile_){
		this(new RTextInputStream(fname1, fname2, maxReads), outname_, minChrom, maxChrom, maxReads, sitesfile_);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
	}
	
	public GenerateVarlets2(RTextInputStream stream_, String outname_, byte minChrom, byte maxChrom, long maxReads, String sitesfile_){
		sitesfile=sitesfile_;
		stream=stream_;
		outname=outname_;
		assert(outname==null || outname.contains("#")) : "Output file name must contain the character '#' to be used for key number.";
		makeKeyMap();
		
		cris=(USE_CRIS ? new ConcurrentLegacyReadInputStream(stream, maxReads) : null);
		if(CONDENSE_SNPS){assert(!SPLIT_SUBS);}
	}
	
	public void finish(){
		
		ArrayList<Long> keys=new ArrayList<Long>();
		keys.addAll(keymap.keySet());
		Shared.sort(keys);
		for(long k : keys){
			ArrayList<Varlet> vars=keymap.remove(k);
			if(!vars.isEmpty()){writeList(vars);}
		}
		
		if(cris!=null){ReadWrite.closeStream(cris);}
		else{stream.close();}
		
	}
	
	public void process(){
		
		Timer t=new Timer();
		
		if(sitesfile!=null){
			sitemap=loadSites(sitesfile);
		}
		
		cris.start();
		ProcessThread[] threadHandles=new ProcessThread[THREADS];
		for(int i=0; i<THREADS; i++){
			threadHandles[i]=new ProcessThread();
			threadHandles[i].start();
		}
		
		long varsMade=0;
		long norefsMade=0;
		long snpMade=0;
		long delMade=0;
		long subnMade=0;
		long subdMade=0;
		long subiMade=0;
		long insMade=0;
		long deltaLen=0;
		long sitesProcessed=0;
		long readsProcessed=0;
		
		for(int i=0; i<threadHandles.length; i++){
			ProcessThread pt=threadHandles[i];
			while(!pt.finished()){
				synchronized(pt){
					try {
						pt.wait(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			varsMade+=pt.varsMade;
			norefsMade+=pt.norefsMade;
			snpMade+=pt.snpMade;
			delMade+=pt.delMade;
			subnMade+=pt.subnMade;
			subdMade+=pt.subdMade;
			subiMade+=pt.subiMade;
			insMade+=pt.insMade;
			deltaLen+=pt.deltaLen;
			sitesProcessed+=pt.sitesProcessed;
			readsProcessed+=pt.readsProcessed;
		}
		
		finish();
		
		t.stop();

		System.out.println("\nOutput variations count");
		System.out.println("Total (minus no-ref):  \t"+(varsMade-norefsMade));
		System.out.println("Deletions:             \t"+(delMade));
		System.out.println("D-type subs:           \t"+(subdMade));
		System.out.println("Insertions:            \t"+(insMade));
		System.out.println("I-type subs:           \t"+(subiMade));
		System.out.println("Snps:                  \t"+(snpMade));
		System.out.println("N-type subs:           \t"+(subnMade));
		System.out.println("No-refs:               \t"+(norefsMade));
		System.out.println("Delta Length:          \t"+(deltaLen));
		System.out.println("Sites Processed:       \t"+(sitesProcessed));
		System.out.println("Reads Processed:       \t"+(readsProcessed));
		System.out.println();
		System.out.println("Time:\t"+t);
	}
	
	
	/**
	 * @param sitesfile2
	 * @return
	 */
	private static final HashMap<Long, SiteR> loadSites(String fname) {
		HashMap<Long, SiteR> map=new HashMap<Long, SiteR>(4096);
		TextFile tf=new TextFile(fname, false);
		
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			SiteScoreR[] array=CalcCoverageFromSites.toSites(s);
			for(SiteScoreR ssr : array){
				SiteR sr=new SiteR(ssr);
				Long key=sr.idPairnum;
				
				SiteR head=map.get(key);
				sr.next=head;
				map.put(key, sr);
			}
			
		}
		return map;
	}
	

	private void writeList(ArrayList<Varlet> list){
		assert(list!=null && list.size()>0);
		long key=key(list.get(0).chromosome, list.get(0).beginLoc);
		String fname=fname(key, outname);
		boolean allowSubprocess=false;
		OutputStream os=ReadWrite.getOutputStream(fname, true, true, allowSubprocess);
		PrintWriter pw=new PrintWriter(os);
		
		
		for(Varlet v : list){
			pw.println(v.toText());
		}
		ReadWrite.finishWriting(pw, os, fname, allowSubprocess);
	}
	 
	
	private final class ProcessThread extends Thread {
		
		public ProcessThread(){
		}
		
		@Override
		public void run(){
			
			final boolean processReads=true;
			if(!processReads){System.err.println("Warning: Skipping read processing.");}
			
			if(cris!=null){
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				while(!terminate && reads!=null && reads.size()>0){
					if(processReads){processReads(reads);}
					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				cris.returnList(ln);
			}else{
				ArrayList<Read> reads=stream.nextList();
				while(!terminate && reads!=null && reads.size()>0){
					if(processReads){processReads(reads);}
					reads=stream.nextList();
				}
			}
			
			finished=true;
			synchronized(this){this.notifyAll();}
		}
		
		private void processReads(ArrayList<Read> reads){

			if(sitemap==null){
				for(Read r : reads){
					Read r2=r.mate;
					assert(r2==null || r.mate.mate==r);

					if(r2==null){
						processRead(r);
					}else{
						if(!TOSS_SOLO1 || r.paired()){processRead(r);}
						if(!TOSS_SOLO2 || r2.paired()){processRead(r2);}
					}
				}
			}else{
				for(Read r : reads){
					Read r2=r.mate;
					assert(r2==null || r.mate.mate==r);
					
					if(r2==null){
						multiprocessRead(r);
					}else{
						if(!TOSS_SOLO1 || r.paired()){multiprocessRead(r);}
						if(!TOSS_SOLO2 || r2.paired()){multiprocessRead(r2);}
					}
				}
			}
		}
		
		@Deprecated
		private void multiprocessRead_old(Read r){
			long key=r.numericID;
			if((r.pairnum()&1)==1){
				key=-key;
				assert(key<0);
			}
			if(true){throw new RuntimeException("Deprecated.");}
			ArrayList<SiteScoreR> alssr=null;//sitemap.get(key);
			if(alssr==null){return;}

			
			for(SiteScoreR ssr : alssr){
				SiteScore ss=find(ssr, r.sites);
				assert(ss!=null) : "\nCan't find ssr "+ssr+" in read\n"+r+"\n";
				
				r.clearSite();
				r.setFromSite(ss);
				r.match=null;
				
				r.setPaired(ss.pairedScore>0);
				r.setPerfect(ss.perfect);
				r.setRescued(ss.rescued);
				
				processRead(r);
			}
		}
		
		private void multiprocessRead(Read r){
			long key=r.numericID;
			if((r.pairnum()&1)==1){
				key=-key;
				assert(key<0);
			}
			
			
			SiteR head=sitemap.get(key);
			if(head!=null){readsProcessed++;}

//			assert(head==null) : "\n"+r.pairnum()+", "+key+",\n"+r.list+",\n"+r.mate.list+"\n"+head.toTextRecursive(null)+"\n";
			
			while(head!=null){
				SiteScore ss=find(head, r.sites);
				assert(ss!=null) : "\nCan't find sr "+head+" in read\n"+r+"\n";
				
				r.clearSite();
				r.setFromSite(ss);
				r.match=null;
				
				r.setPaired(ss.pairedScore>0);
				r.setPerfect(ss.perfect);
				r.setRescued(ss.rescued);
				
				processRead(r);
				SiteR old=head;
				head=old.next;
				old.next=null; //Clears up memory.
			}
		}
		
		/**
		 * @param ssr
		 * @param list
		 * @return
		 */
		private SiteScore find(SiteScoreR ssr, ArrayList<SiteScore> list) {
			for(SiteScore ss : list){
				if(ssr.equals(ss)){return ss;}
			}
			return null;
		}
		
		private SiteScore find(SiteR sr, ArrayList<SiteScore> list) {
			for(SiteScore ss : list){
				if(sr.equals(ss)){return ss;}
			}
			return null;
		}
		

		private void processRead(Read r){
			sitesProcessed++;
			
			boolean flag=false;
			if(false && (/*r.numericID==30719442 ||  r.numericID==107055007  || */ r.numericID==42829556) /*&& r.length()<=35*/){
				System.err.println("Processing read:");
				System.err.println("\n"+r.toText(false));
				System.err.println("\n"+r.strand());
				System.err.println("\n");
				System.err.println(new String(r.bases));
				System.err.println(r.match==null ? "null" : new String(r.match));
				System.err.println("\n");
				tcr.verbose=true;
				flag=true;
				System.err.println("Mapped Length: "+(r.stop-r.start+1));
			}
			
			
//			if(r.chrom<1 && r.list!=null && r.list.size()>0){
//				SiteScore ss=r.list.get(0); //Should not be necessary
//				r.start=ss.start;
//				r.stop=ss.stop;
//				r.chrom=ss.chrom;
//				r.setStrand(ss.strand);
//			}
			assert((r.chrom>=1)==r.mapped()) : r.toText(false);
			if(!r.mapped()){//Unmapped.
				assert(r.sites==null || r.sites.isEmpty()) : r.toText(false);
				return;
			}
			if(r.invalid()){return;} //Probably trimmed too short to be used.
			
			if(r.match!=null){
				if(r.perfect()){//Hopefully this will be set correctly...
					assert(TranslateColorspaceRead.perfectMatch(r.match));
					return;
				}else if(TranslateColorspaceRead.perfectMatch(r.match)){
					return;
				}
			}
			
			if(flag){
				System.err.println("r.match = "+(r.match==null ? null : new String(r.match)));
				System.err.println("Mapped Length: "+(r.stop-r.start+1));
			}
//			if(r.match!=null){
//				for(int i=0; i<r.match.length; i++){
//					if(r.match[i]=='I'){
//						r.match=null;
//						if(flag){System.err.println("nullified match string");}
//						break;
//					}
//				}
//			}
			
//			r.match=null; //TODO - why are some match strings backwards?
			if(r.match==null){
				if(flag){
					System.err.println("realigning match string");
					System.err.println("Mapped Length: "+(r.stop-r.start+1));
				}
				tcr.realign_new(r, 20, true, 0, false); //Also generates the match string
				if(TranslateColorspaceRead.perfectMatch(r.match)){return;}
				if(flag){
					System.err.println("new match string:\n"+(r.match==null ? null : new String(r.match)));
					System.err.println("Mapped Length: "+(r.stop-r.start+1));
				}
			}
			r.errors=r.estimateErrors();
			
			if(r.match==null){
				System.err.println("Could not align read "+r.numericID);
				return;
			}else if(r.match[0]=='X'){
				System.err.println("Could not align read "+r.numericID+": "+new String(r.match));
				return;
			}
			
			assert(CONDENSE);
			ArrayList<Varlet> vars=tcr.toVars(r, CONDENSE, CONDENSE_SNPS, SPLIT_SUBS);
			
			if(vars==null){return;}
			
//			if(r.numericID==36858949){
//				System.err.println(r.toText(false));
//				System.err.println(r.copies);
//				System.err.println(r.mate.toText(false));
//				System.err.println(r.mate.copies);
//				System.err.println();
//
//				for(Varlet v : vars){
//					System.err.println(v.toText());
//					System.err.println(v.numReads);
//				}
//				assert(false);
//			}
			
			for(Varlet v : vars){
				if(v.endDist>=MIN_END_DIST){
					assert(v.numUniqueReads==1);
					assert(v.numSemiUniqueReads==1);
					assert(v.numPlusReads1+v.numMinusReads1+v.numPlusReads2+v.numMinusReads2==1);
					assert(v.numReads>=1);
					//				assert(!TranslateColorspaceReadPacBio.COUNT_DUPLICATES_WHEN_MAKING_VARLETS || v.numReads==1);
					assert(v.numReads==r.copies);
					assert(v.readLen==r.length());
					varsMade++;
					if(v.varType==Variation.NOREF){norefsMade++;}
					else if(v.varType==Variation.SNP){snpMade++;}
					else if(v.varType==Variation.DEL){delMade++;}
					else if(v.varType==Variation.INS){insMade++;}
					else if(v.varType==Variation.DELINS){
						int a=v.lengthRef();
						int b=v.lengthVar();
						if(a==b){subnMade++;}
						else if(a>b){subdMade++;}
						else{subiMade++;}
					}
					deltaLen+=v.lengthDif();
					addVar(v);
				}
			}
//			System.out.println(varsMade+", "+norefsMade);
		}
		
		/** TODO: Synchronize once per read, not once per varlet */
		private void addVar(Varlet v){
			long key=key(v.chromosome, v.beginLoc);
			ArrayList<Varlet> list=keymap.get(key);
			assert(list!=null) : "\nCan't find "+key+" in "+keymap.keySet()+"\n";
			synchronized(list){
				list.add(v);
				if(list.size()>=WRITE_BUFFER){

					if(MERGE_EQUAL_VARLETS){
						mergeEqualVarlets(list);
					}else{
						Shared.sort(list);
					}

					writeList(list);
					list.clear();
				}
			}
		}
		
		private void mergeEqualVarlets(ArrayList<Varlet> vars){
			
			Shared.sort(vars);
			ArrayList<Varlet> list=new ArrayList<Varlet>(8);
			for(int i=0; i<vars.size(); i++){
				Varlet a=vars.get(i);
				vars.set(i, null);
				Varlet b=(list.isEmpty() ? null : list.get(0));
				if(b==null || a.equals(b)){
					list.add(a);
				}else{//purge
					Varlet c=StackVariations.mergeEqualVarlets(list);
					vars.set(i-1, c);
					list.clear();
					list.add(a);
				}
			}
			if(!list.isEmpty()){
				Varlet c=StackVariations.mergeEqualVarlets(list);
				vars.set(list.size()-1, c);
			}
			Tools.condenseStrict(vars);
		}

		protected boolean finished(){return finished;}
		protected void terminate(){terminate=true;}
		
		private final TranslateColorspaceRead tcr=new TranslateColorspaceRead(PAC_BIO_MODE ?
				new MultiStateAligner9ts(ALIGN_ROWS, ALIGN_COLUMNS) :  new MultiStateAligner9ts(ALIGN_ROWS, ALIGN_COLUMNS));
		private boolean finished=false;
		private boolean terminate=false;
		private long varsMade=0;
		private long norefsMade=0;
		private long snpMade=0;
		private long delMade=0;
		private long subnMade=0;
		private long subdMade=0;
		private long subiMade=0;
		private long insMade=0;
		private long deltaLen=0;
		private long sitesProcessed=0;
		private long readsProcessed=0;
		
		
	}
	

	protected static final long key(int chrom, int start){
		long k=((long)chrom<<32)+(Tools.max(start, 0))/BLOCKSIZE;
		return k;
	}
	
	protected static final long[] keys(final int chrom){
		int lim=(Data.chromLengths[chrom]+1000)/BLOCKSIZE;
		long[] keys=new long[lim+1];
		for(int i=0; i<=lim; i++){
			long key=key(chrom, i*BLOCKSIZE);
			keys[i]=key;
		}
		return keys;
	}
	
	protected static final String fname(long key, String outname){
		if(outname==null){outname="GV2TempFile_#.txt";}
		assert(outname.contains("#")) : outname;
		assert(!outname.endsWith(".gz") && !outname.endsWith(".zip") && !outname.endsWith(".bz2")) : outname;
		return outname.replace("#", "b"+Data.GENOME_BUILD+"_"+key);
	}
	
	private final void makeKeyMap(){
		final String header=Varlet.textHeader()+"\n";
		keymap=new HashMap<Long, ArrayList<Varlet>>();
		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			long[] keys=keys(chrom);
			for(long key : keys){
				keymap.put(key, new ArrayList<Varlet>(WRITE_BUFFER));
				ReadWrite.writeString(header, fname(key, outname), false);
			}
		}
	}
	
	private HashMap<Long, ArrayList<Varlet>> keymap;

	public final String outname;
	public final String sitesfile;
//	private HashMap<Long, ArrayList<SiteScoreR>> sitemap=null;
	private HashMap<Long, SiteR> sitemap=null;
	private final RTextInputStream stream;
	private final ConcurrentLegacyReadInputStream cris;
	
	
	public static boolean USE_CRIS=true; //Similar speed either way.  "true" may be better with many threads.
	
	public static int THREADS=7;
	public static int WRITE_BUFFER=200000; //Bigger number uses more memory, for less frequent writes.

	public static boolean CONDENSE=true;
	public static boolean CONDENSE_SNPS=true;
	public static boolean SPLIT_SUBS=false;

	public static boolean TOSS_SOLO1=false;
	public static boolean TOSS_SOLO2=false;
	
	public static boolean MERGE_EQUAL_VARLETS=false;
	public static boolean PAC_BIO_MODE=true;
	public static int ALIGN_ROWS=2020;
	public static int ALIGN_COLUMNS=3000;
	
	public static long MAX_READS=-1;
	public static final int MIN_END_DIST=4;
	public static int BLOCKSIZE=1000000;
	
}
