package pacbio;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import dna.Data;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.SiteScoreR;

/**
 * @author Brian Bushnell
 * @date Aug 2, 2012
 *
 */
public class SortSites {
	
	
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
	
	public static void stack(String fname1, String outname, String tempname){

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
	
	protected static final long key(SiteScoreR ssr){
		return (POSITIONMODE ? poskey(ssr.chrom, ssr.start) : idkey(ssr.numericID));
	}
	
	protected static final long poskey(int chrom, int start){
		long k=((long)chrom<<32)+(Tools.max(start, 0))/BLOCKSIZE;
		return k;
	}
	
	protected static final long idkey(long id){
		long k=id/BLOCKSIZE;
		return k;
	}
	
	protected static final String fname(long key, String outname){
		if(outname==null){outname=DEFAULT_TEMP_PATTERN;}
		assert(outname.contains("#")) : outname;
		return outname.replace("#", "b"+Data.GENOME_BUILD+"_"+key);
	}
	
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

	private static final void finishByPosition(TextStreamWriter out, ArrayList<Long> keys){
		
		
		
		int chrom=0;
		int loc=INTERVAL;
		String tab="";
		StringBuilder sb=new StringBuilder(4000);
		
		for(Long k : keys){
			TextStreamWriter tsw=wmap.get(k);
			String fname=fname(k, tempname);
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
	
	private static final void finishByID(TextStreamWriter out, ArrayList<Long> keys){
		
		long id=0;
		int pairnum=0;
		String tab="";
		StringBuilder sb=new StringBuilder(4000);
		
		for(Long k : keys){
			TextStreamWriter tsw=wmap.get(k);
			String fname=fname(k, tempname);
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
	
	public static int INTERVAL=200;
	public static int BLOCKSIZE=8000000;
	public static long sitesRead=0;
	public static long sitesWritten=0;
	public static long perfectWritten=0;
	public static long semiperfectWritten=0;
	public static boolean DELETE_TEMP=true;
	public static final String DEFAULT_TEMP_PATTERN="SortSitesByIDTempFile_#.txt.gz";
	public static String tempname=null;
	public static boolean POSITIONMODE=false; //False means sort by ID
	public static boolean IGNORE_PERFECT_SITES=false; //Don't process perfect mappings, since they can't yield varlets.
	
}
