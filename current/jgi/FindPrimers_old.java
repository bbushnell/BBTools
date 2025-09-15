package jgi;

import java.io.File;
import java.util.ArrayList;

import align2.MSA;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2014
 *
 */
public class FindPrimers_old {

	public static void main(String[] args){
		Timer t=new Timer();
		FindPrimers_old x=new FindPrimers_old(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public FindPrimers_old(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		float cutoff_=0;
		String literal_=null;
		String ref_=null;
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("replicate") || a.equals("expand")){
				replicateAmbiguous=Parse.parseBoolean(b);
			}else if(a.equals("literal")){
				literal_=b;
			}else if(a.equals("cutoff") || a.equals("minid")){
				cutoff_=Float.parseFloat(b);
				if(cutoff_>1){cutoff_=cutoff_/100;}
				assert(cutoff_>=0 && cutoff_<=1) : "Cutoff should range from 0 to 1";
			}else if(a.equals("primer") || a.equals("query") || a.equals("ref")){
				if(new File(b).exists()){ref_=b;}
				else{literal_=b;}
			}else if(a.equals("msa")){
				msaType=b;
			}else if(a.equals("columns")){
				columns=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
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
		cutoff=cutoff_;
		
		if(ref_!=null){
			ArrayList<Read> list=FastaReadInputStream.toReads(ref_, FileFormat.FASTA, -1);
			int max=0;
			queries=new ArrayList<Read>();
			for(int i=0; i<list.size(); i++){
				Read r=list.get(i);
				max=Tools.max(max, r.length());
				queries.add(r);
				if(rcomp){
					r=r.copy();
					r.reverseComplement();
					r.id="r_"+r.id;
					r.setStrand(1);
					queries.add(r);
				}
			}
			maxqlen=max;
		}else if(literal_!=null){
			int max=0;
			String[] s2=literal_.split(",");
			queries=new ArrayList<Read>();
			for(int i=0; i<s2.length; i++){
				Read r=new Read(s2[i].getBytes(), null, "query", i);
				max=Tools.max(max, r.length());
				queries.add(r);
				
			}
			maxqlen=max;
		}else{
			queries=null;
			maxqlen=0;
		}
		
		if(replicateAmbiguous){
			queries=Tools.replicateAmbiguous(queries, 1);
		}
		if(rcomp){
			for(int i=0, max=queries.size(); i<max; i++){
				Read r=queries.get(i).copy();
				r.reverseComplement();
				r.id="rquery";
				r.setStrand(1);
				queries.add(r);
			}
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		
		final ByteStreamWriter bsw;
		if(out1!=null){

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		
		MSA msa=MSA.makeMSA(maxqlen+3, columns, msaType);
		
		long readsProcessed=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					
					if(r.length()+2>msa.maxColumns){
						msa=MSA.makeMSA(maxqlen+3, r.length()+2+r.length()/2, "MultiStateAligner11ts");
					}
					final int a=0, b=r.length()-1;
					int[] max;
					
					SiteScore bestSite=null;
					Read bestQuery=null;
					for(int qnum=0; qnum<queries.size(); qnum++){
						final Read query=queries.get(qnum);
						max=msa.fillLimited(query.bases, r.bases, a, b, -9999, null);
						if(max!=null){
							int[] score=msa.score(query.bases, r.bases, a, b, max[0], max[1], max[2], false);
							SiteScore ss=new SiteScore(1, query.strand(), score[1], score[2], 1, score[0]);
							if(bestSite==null || ss.quickScore>bestSite.quickScore){
								bestQuery=query;
								ss.setSlowScore(ss.quickScore);
								ss.score=ss.quickScore;
								ss.match=msa.traceback(query.bases, r.bases, a, b, score[3], score[4], score[5], false);
								bestSite=ss;
							}
						}
					}
					
					if(bsw!=null && bestSite!=null){
						bsw.println(toBytes(r, bestQuery, bestSite));
					}
					
					readsProcessed++;
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		if(bsw!=null){bsw.poisonAndWait();}
		ReadWrite.closeStreams(cris);
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}
	

	/*--------------------------------------------------------------*/
	
	private ByteBuilder toBytes(Read r, Read query, SiteScore ss){
		ByteBuilder bb=new ByteBuilder(80);
		float f=Read.identity(ss.match);
		if(f<cutoff){return bb;}
		bb.append(query.id).append('\t');
		bb.append(makeFlag(ss)).append('\t');
		bb.append(r.id.replace('\t', '_')).append('\t');
		bb.append(ss.start+1).append('\t');
		bb.append(Tools.max(ss.score/query.length(), 4)).append('\t');
		String cigar=SamLine.toCigar14(ss.match, ss.start, ss.stop, r.length(), query.bases);
		if(cigar==null){bb.append('*').append('\t');}else{bb.append(cigar).append('\t');}
		bb.append('0').append('\t');
		bb.append('*').append('\t');
		bb.append('0').append('\t');
		
		bb.append(query.bases).append('\t');
		bb.append('*').append('\t');
		
		bb.append(Tools.format("YI:f:%.2f", (100*f)));
		
		return bb;
	}
	
	public static int makeFlag(SiteScore ss){
		int flag=0;
		if(ss.strand()==Shared.MINUS){flag|=0x10;}
		return flag;
	}
	
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private final float cutoff;
	private boolean rcomp=true;
	private boolean replicateAmbiguous=true;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private ArrayList<Read> queries;
	private final int maxqlen;
	private int columns=2000;
	private String msaType="MultiStateAligner11ts";
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
