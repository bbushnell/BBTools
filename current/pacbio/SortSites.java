package pacbio;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import dna.Data;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.SiteScoreR;

/**
 * External-sorts tab-delimited SiteScoreR lists either by genomic position or
 * by read numeric ID/pair number, optionally dropping perfect sites before the
 * temporary-file split.
 * @author Brian Bushnell
 * @date Aug 2, 2012
 */
public class SortSites{
	
	/**
	 * Parses mode, block and temporary-file flags, sorts the input site table, and
	 * prints global site counters to standard output.
	 * @param args Input path, output path, and optional flag=value settings
	 */
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		String tempname=null;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("genome") || a.equals("build")){
				Data.setGenome(Integer.parseInt(b)); //Not needed
			}else if(a.equals("tempname")){
				tempname=b;
			}else if(a.equals("deletefiles") || a.startsWith("deletetemp") || a.equals("delete")){
				DELETE_TEMP=(Parse.parseBoolean(b));
			}else if(a.equals("mode")){
				POSITIONMODE=(b.contains("position") || b.contains("location"));
			}else if(a.equals("blocksize")){
				BLOCKSIZE=(Integer.parseInt(b));
			}else if(a.equals("ignoreperfect")){
				IGNORE_PERFECT_SITES=(Parse.parseBoolean(b));
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		if(POSITIONMODE){
			System.out.println("Sorting by position.");
		}else{
			System.out.println("Sorting by ID.");
		}
		
		stack(args[0], args[1], tempname);
		assert(sitesRead==sitesWritten || (sitesRead>=sitesWritten && IGNORE_PERFECT_SITES));
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	/**
	 * Splits parsed SiteScoreR records into temporary buckets, then finishes the
	 * selected global order into one output file. Perfect sites are counted as
	 * read but skipped when {@link #IGNORE_PERFECT_SITES} is set.
	 * @param fname1 Input tab-delimited SiteScoreR list
	 * @param outname Output sorted site-list path
	 * @param tempname Optional temporary filename pattern copied to the static
	 * field used by bucket writers
	 */
	public static void stack(String fname1, String outname, String tempname){

		SortSites.tempname=tempname;

		TextFile tf=new TextFile(fname1, false);

		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){

			SiteScoreR[] array=SiteScoreR.fromTextArray(s);
			sitesRead+=array.length;

			for(SiteScoreR ssr : array){
				if(!ssr.perfect || !IGNORE_PERFECT_SITES){
					write(ssr);
				}
			}
		}
		tf.close();

		System.out.println("Finished reading");
		System.out.println("Read "+sitesRead+" sites.");

		finish(outname);
		System.out.println("Wrote "+sitesWritten+" sites.");
		System.out.println("Wrote "+perfectWritten+" perfect sites.");
		System.out.println("Wrote "+semiperfectWritten+" semiperfect sites.");
		wmap.clear();
	}
	
	/** Queues one site to its temp-file bucket, opening the bucket writer on demand. */
	private static void write(SiteScoreR ssr){
		long key=key(ssr);
		TextStreamWriter tsw=wmap.get(key);
		if(tsw==null){
			String fname=fname(key, tempname);
			tsw=new TextStreamWriter(fname, true, false, false);
			tsw.start();
			wmap.put(key, tsw);
		}
		tsw.print(ssr.toText().append('\n'));
	}
	
	/** Returns the temp-file bucket key for position or ID sorting mode. */
	protected static final long key(SiteScoreR ssr){
		return (POSITIONMODE ? poskey(ssr.chrom, ssr.start) : idkey(ssr.numericID));
	}
	
	/** Returns the genomic-position bucket for one chromosome and start offset. */
	protected static final long poskey(int chrom, int start){
		long k=((long)chrom<<32)+(Tools.max(start, 0))/BLOCKSIZE;
		return k;
	}
	
	/** Returns the read-ID bucket for one numeric read ID. */
	protected static final long idkey(long id){
		long k=id/BLOCKSIZE;
		return k;
	}
	
	/** Expands a temp-name pattern with the current genome build and bucket key. */
	protected static final String fname(long key, String outname){
		if(outname==null){outname=DEFAULT_TEMP_PATTERN;}
		assert(outname.contains("#")) : outname;
		return outname.replace("#", "b"+Data.GENOME_BUILD+"_"+key);
	}
	
	//n Historical traced-clean note from a prior pass: [SortSites] TRACED CLEAN (no .sh, no callers — invoked by MakePacBioScript scripts). External-sort-via-temp-files, same
	//n architecture as StackSites2. Correct: temp writers poisoned (finish L140) THEN joined per-block (finishBy* L164/L232)
	//n before their files are read — proper flush ordering. finishByID grouping (`numericID>id||pairnum>pairnum` + else-assert
	//n `id==id&&pairnum==pairnum`) is SOUND: IDCOMP sorts (id,pairnum) monotonically and a given id never spans blocks
	//n (idkey=id/BLOCKSIZE), so within-id pairnum is non-decreasing → the else-assert can't fire wrongly. poskey/idkey packing
	//n correct. Only nits: args[0..1] unguarded (AIOOBE<2); the unconditional sb.append('\n')+print at finishBy* end emits a
	//n lone newline when the last block was empty (cosmetic). sitesRead==sitesWritten sanity assert (main L68) is consistent.
	/**
	 * Poisons all temp writers, then merges temp buckets in sorted bucket-key order.
	 * The final output writer is joined, but its terminal error state is ignored.
	 * @param outname Output sorted site-list path
	 */
	private static final void finish(String outname){
		TextStreamWriter out=new TextStreamWriter(outname, true, false, false);
		out.start();
		ArrayList<Long> keys=new ArrayList<Long>(wmap.size());
		keys.addAll(wmap.keySet());
		Shared.sort(keys);
		for(Long k : keys){
			TextStreamWriter tsw=wmap.get(k);
			tsw.poison();
		}
		
		if(POSITIONMODE){
			finishByPosition(out, keys);
		}else{
			finishByID(out, keys);
		}

		out.poisonAndWait();
	}

	/**
	 * Joins each temp writer, reads and sorts each bucket by genomic coordinate,
	 * and emits tab-separated records grouped into {@link #INTERVAL}-base spans.
	 * @param out Output writer
	 * @param keys Sorted temp bucket keys
	 */
	private static final void finishByPosition(TextStreamWriter out, ArrayList<Long> keys){
		int chrom=0;
		int loc=INTERVAL;
		String tab="";
		StringBuilder sb=new StringBuilder(4000);
		
		for(Long k : keys){
			TextStreamWriter tsw=wmap.get(k);
			String fname=fname(k, tempname);
			for(int i=0; i<50 && tsw.isAlive(); i++){
				try{
					tsw.join(20000);
				}catch(InterruptedException e){
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
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){list.add(SiteScoreR.fromText(s));}
			tf.close();
			if(DELETE_TEMP){
				new File(fname).delete();
			}
			Shared.sort(list, SiteScoreR.PCOMP);
			
			final int lim=list.size();
			for(int i=0; i<lim; i++){
				SiteScoreR ssr=list.get(i);
				sitesWritten++;
				if(ssr.semiperfect){semiperfectWritten++;}
				if(ssr.perfect){perfectWritten++;}
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
	}
	
	/**
	 * Joins each temp writer, reads and sorts each bucket by numeric ID/pairnum,
	 * and emits tab-separated records grouped by read mate.
	 * @param out Output writer
	 * @param keys Sorted temp bucket keys
	 */
	private static final void finishByID(TextStreamWriter out, ArrayList<Long> keys){
		long id=0;
		int pairnum=0;
		String tab="";
		StringBuilder sb=new StringBuilder(4000);
		
		for(Long k : keys){
			TextStreamWriter tsw=wmap.get(k);
			String fname=fname(k, tempname);
			for(int i=0; i<50 && tsw.isAlive(); i++){
				try{
					tsw.join(20000);
				}catch(InterruptedException e){
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
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){list.add(SiteScoreR.fromText(s));}
			tf.close();
			if(DELETE_TEMP){
				new File(fname).delete();
			}
			Shared.sort(list, SiteScoreR.IDCOMP);
			
			final int lim=list.size();
			for(int i=0; i<lim; i++){
				SiteScoreR ssr=list.get(i);
				sitesWritten++;
				if(ssr.semiperfect){semiperfectWritten++;}
				if(ssr.perfect){perfectWritten++;}
				list.set(i, null);
				if(ssr.numericID>id || ssr.pairnum>pairnum){
					if(sb.length()>0){//Purge to disk
						sb.append('\n');
						out.print(sb.toString());
						sb.setLength(0);
					}
					id=ssr.numericID;
					pairnum=ssr.pairnum;
					tab="";
				}else{
					assert(ssr.numericID==id && ssr.pairnum==pairnum);
				}
				sb.append(tab);
				sb.append(ssr.toText());
				tab="\t";
			}
			
		}
		
		
		sb.append('\n');
		out.print(sb.toString());
	}

	private static final HashMap<Long, TextStreamWriter> wmap=new HashMap<Long, TextStreamWriter>();
	
	/** Genomic span per position-mode output group */
	public static int INTERVAL=200;
	/** Number of bases or read IDs per temp-file bucket */
	public static int BLOCKSIZE=8000000;
	/** Number of parsed SiteScoreR records */
	public static long sitesRead=0;
	/** Number of emitted SiteScoreR records */
	public static long sitesWritten=0;
	/** Number of emitted perfect SiteScoreR records */
	public static long perfectWritten=0;
	/** Number of emitted semiperfect SiteScoreR records */
	public static long semiperfectWritten=0;
	/** Delete temp buckets after each bucket has been read */
	public static boolean DELETE_TEMP=true;
	/** Default temp-file pattern; # is replaced with build and bucket key */
	public static final String DEFAULT_TEMP_PATTERN="SortSitesByIDTempFile_#.txt.gz";
	/** Temp-file pattern supplied on the command line */
	public static String tempname=null;
	/** False sorts by read ID; true sorts by genomic position */
	public static boolean POSITIONMODE=false;
	/** Do not write perfect mappings, since they can't yield varlets */
	public static boolean IGNORE_PERFECT_SITES=false;
	
}
