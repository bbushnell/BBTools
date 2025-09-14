package pacbio;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import align2.MultiStateAligner9PacBio;
import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentLegacyReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;
import stream.SiteScoreR;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Jul 16, 2012
 *
 */
public class StackSites2 {
	
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}

		Timer t=new Timer();
		
		String tempname=null;
		Data.GENOME_BUILD=-1;
		
		for(int i=4; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("genome") || a.equals("build")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("tempname")){
				tempname=b;
			}else if(a.equals("deletefiles") || a.startsWith("deletetemp") || a.equals("delete")){
				DELETE_TEMP=(Parse.parseBoolean(b));
			}else if(a.equals("blocksize")){
				BLOCKSIZE=(Integer.parseInt(b));
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(Data.GENOME_BUILD<0){throw new RuntimeException("Please specify genome build.");}
		
		stack(args[0], args[1], args[2], args[3], tempname);
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	public static void stack(String fname1, String fname2, String outname, String pcovoutname, String tempname){
		assert(pcovoutname.contains("#"));
		final RTextInputStream rtis=new RTextInputStream(fname1, (fname2==null || fname2.equals("null") ? null : fname2), -1);
		final ConcurrentLegacyReadInputStream cris=new ConcurrentLegacyReadInputStream(rtis, -1);
		
		cris.start();
		System.err.println("Started cris");
		final boolean paired=cris.paired();
		System.err.println("Paired: "+paired);
		
		final ArrayList<CoverageArray> pcov;
		final ArrayList<CoverageArray> truePcov;
		final ArrayList<CoverageArray> cov;
		
		{
			int len=(Data.GENOME_BUILD<0 ? 8 : Data.numChroms+1);
			
			pcov=new ArrayList<CoverageArray>(len);
			truePcov=new ArrayList<CoverageArray>(len);
			cov=new ArrayList<CoverageArray>(len);
			
			System.out.println("len="+len+"; Data.numChroms="+Data.numChroms);
			
			pcov.add(null);
			truePcov.add(null);
			cov.add(null);
			
			for(int i=1; i<len; i++){
				if(Data.GENOME_BUILD<0){
					pcov.add(new CoverageArray2(-1, 500));
					truePcov.add(new CoverageArray2(-1, 500));
					cov.add(new CoverageArray2(-1, 500));
				}else{
					pcov.add(new CoverageArray2(-1, Data.chromLengths[i]+1));
					truePcov.add(new CoverageArray2(-1, Data.chromLengths[i]+1));
					cov.add(new CoverageArray2(-1, Data.chromLengths[i]+1));
				}
			}
		}

		
		final Glob g=new Glob(tempname);
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(paired==(r.mate!=null));
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					readsProcessed++;
					
//					System.out.println("Processing read "+r.numericID);
					
					if(r!=null){
						if(r.sites!=null){
//							System.out.println("Adding "+r.list.size()+" sites.");
							SiteScore original=r.originalSite;
							for(SiteScore ss : r.sites){
								sitesProcessed++;
								
								//TODO: Process perfect coverage
								{
									boolean b=false;
									if(ss.semiperfect){
										b=true;
									}else{//Check for no-refs
										int len=ss.stop-ss.start+1;
										if(len==r.length() && ss.slowScore>=0.5f*MultiStateAligner9PacBio.POINTS_MATCH2){
											b=checkPerfection(ss.start, ss.stop, r.bases, Data.getChromosome(ss.chrom), ss.strand==Shared.MINUS, 0.5f);
										}
									}
									if(b){
										while(pcov.size()<=ss.chrom){
											pcov.add(new CoverageArray2(pcov.size(), Data.chromLengths[pcov.size()]));
											truePcov.add(new CoverageArray2(truePcov.size(), Data.chromLengths[truePcov.size()]));
										}
										CoverageArray ca=pcov.get(ss.chrom);
										CoverageArray tca=truePcov.get(ss.chrom);
										for(int i=ss.start+PCOV_TIP_DIST; i<=ss.stop-PCOV_TIP_DIST; i++){
											ca.increment(i);
										}
										if(ss.perfect){
											for(int i=ss.start; i<=ss.stop; i++){
												tca.increment(i);
											}
										}
									}
									{
										while(cov.size()<=ss.chrom){
											cov.add(new CoverageArray2(cov.size(), Data.chromLengths[cov.size()]));
										}
										CoverageArray ca=cov.get(ss.chrom);
										for(int i=ss.start; i<=ss.stop; i++){
											ca.increment(i);
										}
									}
								}
								
								SiteScoreR ssr=new SiteScoreR(ss, r.length(), r.numericID, (byte)r.pairnum());
								
								if(original!=null){
									ssr.correct=isCorrectHitLoose(ss, original.chrom, original.strand, original.start, original.stop, 40, false);
								}
								
								g.write(ssr);
							}
//							System.out.println(sitesProcessed);
						}
					}
					
					if(r.mate!=null){
						Read r2=r.mate;
						if(r2.sites!=null){
							
							SiteScore original=r2.originalSite;
							for(SiteScore ss : r2.sites){
								sitesProcessed++;
								
								{
									boolean b=false;
									if(ss.semiperfect){
										b=true;
									}else{//Check for no-refs
										int len=ss.stop-ss.start+1;
										if(len==r2.length() && ss.slowScore>=0.5f*MultiStateAligner9PacBio.POINTS_MATCH2){
											b=checkPerfection(ss.start, ss.stop, r2.bases, Data.getChromosome(ss.chrom), ss.strand==Shared.MINUS, 0.5f);
										}
									}
									if(b){
										while(pcov.size()<=ss.chrom){
											pcov.add(new CoverageArray2(pcov.size(), Data.chromLengths[pcov.size()]));
											truePcov.add(new CoverageArray2(truePcov.size(), Data.chromLengths[truePcov.size()]));
										}
										CoverageArray ca=pcov.get(ss.chrom);
										CoverageArray tca=truePcov.get(ss.chrom);
										for(int i=ss.start+PCOV_TIP_DIST; i<=ss.stop-PCOV_TIP_DIST; i++){
											ca.increment(i);
										}
										if(ss.perfect){
											for(int i=ss.start; i<=ss.stop; i++){
												tca.increment(i);
											}
										}
									}
									{
										while(cov.size()<=ss.chrom){
											cov.add(new CoverageArray2(cov.size(), Data.chromLengths[cov.size()]));
										}
										CoverageArray ca=cov.get(ss.chrom);
										for(int i=ss.start; i<=ss.stop; i++){
											ca.increment(i);
										}
									}
								}
								
								SiteScoreR ssr=new SiteScoreR(ss, r2.length(), r2.numericID, (byte)r2.pairnum());
								
								if(original!=null){
									ssr.correct=isCorrectHitLoose(ss, original.chrom, original.strand, original.start, original.stop, 40, false);
								}
								
								g.write(ssr);
							}
						}
					}
					
//					System.out.println(r.toString());
//					assert(r.list!=null);
//					assert(r.list.size()>0);
					
				}
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			System.out.println("Finished reading");
			cris.returnList(ln);
			System.out.println("Returned list");
			ReadWrite.closeStream(cris);
			System.out.println("Closed stream");
			System.out.println("Processed "+readsProcessed+" reads.");
			System.out.println("Processed "+sitesProcessed+" sites.");
		}
		
		
		for(int i=1; i<pcov.size(); i++){
			CoverageArray ca=pcov.get(i);
//			pcov.set(i, null);
			if(ca.maxIndex<.995*ca.arrayLength()){
				ca.resize(ca.maxIndex+1);
			}
			ReadWrite.writeObjectInThread(ca, pcovoutname.replaceFirst("#", ""+i), false);
		}
		
		finish(g, outname, pcov, truePcov, cov);
		System.out.println("Retained  "+sitesOut+" sites.");
		
	}
	
	/** TODO - thread this by chrom */
	private static final void finish(Glob g, String outname, ArrayList<CoverageArray> pcov, ArrayList<CoverageArray> truePcov, ArrayList<CoverageArray> cov){

		
		final TextStreamWriter out=new TextStreamWriter(outname, true, false, false);
		out.start();
		ArrayList<Long> keys=new ArrayList<Long>(g.wmap.size());
		keys.addAll(g.wmap.keySet());
		Shared.sort(keys);
		for(Long k : keys){
			TextStreamWriter tsw=g.wmap.get(k);
			tsw.poison();
		}
		
		
		
		int chrom=0;
		int loc=INTERVAL;
		String tab="";
		StringBuilder sb=new StringBuilder(4000);
		
		for(Long k : keys){
			TextStreamWriter tsw=g.wmap.get(k);
			String fname=Glob.fname(k, g.tempname);
			for(int i=0; i<50 && tsw.isAlive(); i++){
				try {
					tsw.join(20000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(tsw.isAlive()){
					System.err.println("Waiting for tsw "+tsw.fname+" to finish...");
				}
			}
			if(tsw.isAlive()){
				System.err.println(tsw.getClass().getName()+" for "+fname+" refused to die after a long time.");
				assert(false);
			}
			
			TextFile tf=new TextFile(fname, false);
			ArrayList<SiteScoreR> list=new ArrayList<SiteScoreR>(1000);
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				SiteScoreR ssr=SiteScoreR.fromText(s);
				
				assert(pcov.size()>=ssr.chrom) : ssr.chrom+", "+pcov.size()+", "+truePcov.size()+", "+cov.size();
				final int c=ssr.chrom;
				boolean retain=retainSite(ssr, (pcov.size()>c ? pcov.get(c) : FAKE), (truePcov.size()>c ? truePcov.get(c) : FAKE), (cov.size()>c ? cov.get(c) : null));
				if(retain){
					list.add(ssr);
					sitesOut++;
				}
			}
			tf.close();
			if(DELETE_TEMP){
				new File(fname).delete();
			}
			Shared.sort(list, SiteScoreR.PCOMP);
			
			final int lim=list.size();
			for(int i=0; i<lim; i++){
				SiteScoreR ssr=list.get(i);
				list.set(i, null);
				if(ssr.chrom>chrom || ssr.start>=loc){
					if(sb.length()>0){//Purge to disk
						sb.append('\n');
						out.print(sb.toString());
						sb.setLength(0);
					}
					chrom=ssr.chrom;
					loc=ssr.start;
					loc=(loc-(loc%INTERVAL))+INTERVAL;
					assert(loc>ssr.start);
					assert(loc-ssr.start<=INTERVAL);
					assert(loc%INTERVAL==0);
					tab="";
				}
				sb.append(tab);
				sb.append(ssr.toText());
				tab="\t";
			}
			
		}
		
		
		sb.append('\n');
		out.print(sb.toString());
		out.poisonAndWait();
	}
	
	private static boolean retainSite(SiteScoreR ssr, CoverageArray pcov, CoverageArray tpcov, CoverageArray cov){
		if(ssr.semiperfect && !ssr.perfect){return true;} //For tip extension
		assert(cov!=null && cov!=FAKE) : (cov==FAKE)+", "+ssr.chrom;
		
		if(!ssr.semiperfect){ //Typical flawed read
			assert(!ssr.perfect);
			boolean toss=true;
			if(pcov==null || tpcov==null){
				toss=false;
			}else{
				for(int j=ssr.start-PCOV_TIP_DIST; toss && j<=ssr.stop+PCOV_TIP_DIST; j++){
					toss=(pcov.get(j)>=MIN_PCOV_TO_TOSS && tpcov.get(j)>=MIN_PCOV_TO_TOSS);
				}
			}
			if(toss){
				for(int j=ssr.start; j<=ssr.stop; j++){cov.increment(j, -1);}
				return false;
			}
		}

		boolean alwaysLowCov=true;
		boolean alwaysTooPerfect=true;
		boolean onlyPerfect=true;
		
		for(int j=ssr.start; (alwaysLowCov || alwaysTooPerfect || onlyPerfect) && j<=ssr.stop; j++){
			int c=cov.get(j);
			int tp=tpcov.get(j);
			
			alwaysLowCov=alwaysLowCov && c<MIN_COV_TO_RETAIN;
			alwaysTooPerfect=alwaysTooPerfect && c-tp<tp;
			onlyPerfect=onlyPerfect && tp>0;
		}
		
		if(alwaysLowCov || (alwaysTooPerfect && !ssr.semiperfect) || onlyPerfect){
			if(!ssr.semiperfect){
				for(int j=ssr.start; j<=ssr.stop; j++){cov.increment(j, -1);}
			}
			return false;
		}
		
		return true;
	}
	
	private static boolean checkPerfection(int start, int stop, byte[] bases, ChromosomeArray cha, boolean rcomp, float f) {
		
		int noref=0;
		if(rcomp){
			for(int i=0; i<bases.length; i++){
				byte a=AminoAcid.baseToComplementExtended[bases[bases.length-i-1]];
				byte b=cha.get(start+i);
				if(b=='N'){noref++;}
				else if(a!=b){return false;}
			}
		}else{
			for(int i=0; i<bases.length; i++){
				byte a=bases[i];
				byte b=cha.get(start+i);
				if(b=='N'){noref++;}
				else if(a!=b){return false;}
			}
		}
		return bases.length-noref>=f*bases.length;
	}
	
	public static boolean isCorrectHitLoose(SiteScore ss, int trueChrom, byte trueStrand, int trueStart, int trueStop, int thresh, boolean useChrom){
		if((useChrom && ss.chrom!=trueChrom) || ss.strand!=trueStrand){return false;}

		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;
		
		return (Tools.absdif(ss.start, trueStart)<=thresh || Tools.absdif(ss.stop, trueStop)<=thresh);
	}
	
	private static class Glob{
		
		public Glob(String tempPattern_){
			tempname=(tempPattern_ == null ? DEFAULT_TEMP_PATTERN : tempPattern_);
		}
		
		public void write(SiteScoreR ssr){
			long key=key(ssr.chrom, ssr.start);
			TextStreamWriter tsw=wmap.get(key);
			if(tsw==null){
				String fname=fname(key, tempname);
				tsw=new TextStreamWriter(fname, true, false, false);
				tsw.start();
				wmap.put(key, tsw);
			}
			tsw.print(ssr.toText().append('\n'));
		}
		
		protected static final long key(int chrom, int start){
			long k=((long)chrom<<32)+(Tools.max(start, 0))/BLOCKSIZE;
			return k;
		}
		
		protected static final String fname(long key, String outname){
			if(outname==null){outname=DEFAULT_TEMP_PATTERN;}
			assert(outname.contains("#")) : outname;
			return outname.replace("#", "b"+Data.GENOME_BUILD+"_"+key);
		}

		final HashMap<Long, TextStreamWriter> wmap=new HashMap<Long, TextStreamWriter>();
		final String tempname;
		
	}
	
	/** Sites will be written to files, each containing an index range of this size.
	 * Larger means fewer files, but more memory used when reading the files (at a later stage).
	 */
	public static int BLOCKSIZE=8000000;
	
	/** Sites are grouped into intervals (by start location) and treated as an array of arrays.
	 * All sites in an interval are printed as one line of text. */
	public static final int INTERVAL=200;
	public static long readsProcessed=0;
	public static long sitesProcessed=0;
	public static long sitesOut=0;
	public static boolean DELETE_TEMP=true;
	public static final String DEFAULT_TEMP_PATTERN="StackSites2TempFile_#.txt.gz";
	/** Start incrementing coverage this far in from the site tips. */
	public static int PCOV_TIP_DIST=6;

	/** Toss sites from areas with less than this coverage, since they can't be used to call vars */
	public static int MIN_COV_TO_RETAIN=2;
	/** Toss sites from areas with less than this coverage, since they can't be used to call vars */
	public static int MIN_PCOV_TO_TOSS=3;
	
	private static final CoverageArray FAKE=new CoverageArray2(-1, 500);
}
