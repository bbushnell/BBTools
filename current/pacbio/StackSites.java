package pacbio;

import java.util.ArrayList;

import align2.MultiStateAligner9PacBio;
import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
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
public class StackSites {
	
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		for(int i=4; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("genome") || a.equals("build")){
				Data.setGenome(Integer.parseInt(b));
			}
		}
		
		stack(args[0], args[1], args[2], args[3]);
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	public static void stack(String fname1, String fname2, String outname, String pcovoutname){
		assert(pcovoutname.contains("#"));
		RTextInputStream rtis=new RTextInputStream(fname1, (fname2==null || fname2.equals("null") ? null : fname2), -1);
		ConcurrentLegacyReadInputStream cris=new ConcurrentLegacyReadInputStream(rtis, -1);
		
		cris.start();
		System.err.println("Started cris");
		boolean paired=cris.paired();
		System.err.println("Paired: "+paired);
		
		ArrayList<CoverageArray> pcov=new ArrayList<CoverageArray>(8);
		pcov.add(new CoverageArray2(0,1000));
		
		Glob g=new Glob();
		
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
									if(ss.perfect || ss.semiperfect){
										b=true;
									}else{//Check for no-refs
										int len=ss.stop-ss.start+1;
										if(len==r.length() && ss.slowScore>=0.5f*MultiStateAligner9PacBio.POINTS_MATCH2){
											b=checkPerfection(ss.start, ss.stop, r.bases, Data.getChromosome(ss.chrom), ss.strand==Shared.MINUS, 0.5f);
										}
									}
									if(b){
										while(pcov.size()<=ss.chrom){
											pcov.add(new CoverageArray2(pcov.size(), 500));
										}
										CoverageArray ca=pcov.get(ss.chrom);
										for(int i=ss.start; i<=ss.stop; i++){
											ca.increment(i);
										}
									}
								}
								
								SiteScoreR ssr=new SiteScoreR(ss, r.length(), r.numericID, (byte)r.pairnum());
								
								if(original!=null){
									ssr.correct=isCorrectHitLoose(ss, original.chrom, original.strand, original.start, original.stop, 40, false);
								}
								
								g.add(ssr);
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
									if(ss.perfect || ss.semiperfect){
										b=true;
									}else{//Check for no-refs
										int len=ss.stop-ss.start+1;
										if(len==r2.length() && ss.slowScore>=0.5f*MultiStateAligner9PacBio.POINTS_MATCH2){
											b=checkPerfection(ss.start, ss.stop, r2.bases, Data.getChromosome(ss.chrom), ss.strand==Shared.MINUS, 0.5f);
										}
									}
									if(b){
										while(pcov.size()<=ss.chrom){
											pcov.add(new CoverageArray2(pcov.size(), 500));
										}
										CoverageArray ca=pcov.get(ss.chrom);
										for(int i=ss.start; i<=ss.stop; i++){
											ca.increment(i);
										}
									}
								}
								
								SiteScoreR ssr=new SiteScoreR(ss, r2.length(), r2.numericID, (byte)r2.pairnum());
								
								if(original!=null){
									ssr.correct=isCorrectHitLoose(ss, original.chrom, original.strand, original.start, original.stop, 40, false);
								}
								
								g.add(ssr);
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
			System.err.println("Finished reading");
			cris.returnList(ln);
			System.err.println("Returned list");
			ReadWrite.closeStream(cris);
			System.err.println("Closed stream");
			System.err.println("Processed "+readsProcessed+" reads.");
			System.err.println("Processed "+sitesProcessed+" sites.");
		}
		
		
		for(int i=1; i<pcov.size(); i++){
			CoverageArray ca=pcov.get(i);
			pcov.set(i, null);
			ca.resize(ca.maxIndex+1);
			ReadWrite.writeObjectInThread(ca, pcovoutname.replaceFirst("#", ""+i), false);
		}
		

		TextStreamWriter out=new TextStreamWriter(outname, true, false, false);
		out.start();
		
		//This is split by chrom to avoid getting more than 2^31 sites in a single list.
		//Ultimately, it may be better to split output files by chrom and output unsorted, to avoid memory usage.
		for(int i=0; i<g.array.length; i++){
			Shared.sort(g.array[i], SiteScoreR.PCOMP);
			write(g.array[i], out);
			g.array[i]=null;
		}
		
		out.poison();
		
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

	private static void write(ArrayList<SiteScoreR> alsr, TextStreamWriter out){
		if(alsr==null || alsr.size()==0){return;}
		
		int chrom=0;
		int loc=INTERVAL;
		StringBuilder sb=new StringBuilder();
		
		String tab="";
		
		final int lim=alsr.size();
		for(int i=0; i<lim; i++){
			SiteScoreR ssr=alsr.get(i);
			alsr.set(i, null);
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
		
		if(sb.length()>0){//Purge to disk
			sb.append('\n');
			out.print(sb.toString());
			sb.setLength(0);
		}
	}
	
	public static boolean isCorrectHitLoose(SiteScore ss, int trueChrom, byte trueStrand, int trueStart, int trueStop, int thresh, boolean useChrom){
		if((useChrom && ss.chrom!=trueChrom) || ss.strand!=trueStrand){return false;}

		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;
		
		return (Tools.absdif(ss.start, trueStart)<=thresh || Tools.absdif(ss.stop, trueStop)<=thresh);
	}
	
	private static class Glob{
		
		public Glob(){
			array=new ArrayList[8];
			for(int i=0; i<array.length; i++){
				array[i]=new ArrayList<SiteScoreR>();
			}
		}
		
		public void add(SiteScoreR ssr){
			if(ssr.chrom>=array.length){
				int newlen=((int)ssr.chrom*2);
				assert(newlen>array.length);
				ArrayList<SiteScoreR>[] array2=new ArrayList[newlen];
				for(int i=0; i<array.length; i++){array2[i]=array[i];}
				for(int i=array.length; i<array2.length; i++){array2[i]=new ArrayList<SiteScoreR>();}
				array=array2;
			}
			array[ssr.chrom].add(ssr);
		}
		
		public ArrayList<SiteScoreR>[] array;
		
	}
	
	public static final int INTERVAL=200;
	public static long readsProcessed=0;
	public static long sitesProcessed=0;
	
}
