package pacbio;

import java.util.ArrayList;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Timer;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Nov 15, 2012
 *
 */
public class PartitionReads {

	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		boolean verbose=false;
		int ziplevel=-1;
		String in1=null;
		String in2=null;
		long maxReads=-1;
		String outname1=null;
		String outname2=null;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
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
			}else if(a.equals("path") || a.equals("root") || a.equals("tempdir")){
				Data.setPath(b);
			}else if(a.equals("fasta") || a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
				if(b.indexOf('#')>-1){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.startsWith("partition")){
				partitions=Integer.parseInt(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
				System.out.println("Set overwrite to "+overwrite);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("out") || a.equals("out1")){
				if(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none") || split.length==1){
					System.out.println("No output file.");
					outname1=null;
				}else{
					outname1=b;
					assert(!outname1.equalsIgnoreCase(outname2));
				}
			}else if(a.equals("out2")){
				if(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none") || split.length==1){
					outname2=null;
				}else{
					outname2=b;
					assert(!outname2.equalsIgnoreCase(outname1));
				}
			}else if(a.startsWith("verbose")){
				verbose=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		assert(FastaReadInputStream.settingsOK());
		assert(outname1==null || outname1.indexOf('#')>=0 || partitions<2);
		assert(outname2==null || outname2.indexOf('#')>=0 || partitions<2);
		assert(outname1==null || !outname1.equalsIgnoreCase(outname2));
		
		if(in1==null){throw new RuntimeException("Please specify input file.");}
		
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
//			cris.start(); //4567
//			th.start();
		}
		

		TextStreamWriter[] tsw1=new TextStreamWriter[partitions];
		TextStreamWriter[] tsw2=new TextStreamWriter[partitions];
		
		FileFormat ff=FileFormat.testOutput(outname1, FileFormat.FASTQ, null, true, overwrite, append, false);
		fastq=ff.fastq();
		fasta=ff.fasta();
		bread=ff.bread();
		
		for(int i=0; i<partitions; i++){
			tsw1[i]=new TextStreamWriter(outname1.replaceFirst("#", ""+i), overwrite, false, true);
			if(outname2!=null){
				tsw2[i]=new TextStreamWriter(outname2.replaceFirst("#", ""+i), overwrite, false, true);
			}
		}
		
		long reads=process(tsw1, tsw2, cris);
		t.stop();
		System.out.println("Reads: \t"+reads);
		System.out.println("Time:  \t"+t);
	}
	
	public static long process(TextStreamWriter[] tsw1, TextStreamWriter[] tsw2, ConcurrentReadInputStream cris){
		for(TextStreamWriter tsw : tsw1){if(tsw!=null){tsw.start();}}
		for(TextStreamWriter tsw : tsw2){if(tsw!=null){tsw.start();}}
		cris.start(); //4567
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> readlist=ln.list;
		
		final boolean paired=cris.paired();
		
		long x=0;
		final int div=tsw1.length;
		while(!readlist.isEmpty()){
			
			//System.err.println("Got a list of size "+readlist.size());
			for(int i=0; i<readlist.size(); i++){
				Read r=readlist.get(i);
				if(r!=null){
					final Read r2=r.mate;
					final int mod=(int)(x%div);
					
					ByteBuilder a=null, b=null;
					
					if(fastq){
						a=r.toFastq();
						if(paired){b=r2.toFastq();}
					}else if(fasta){
						a=r.toFasta();
						if(paired){b=r2.toFasta();}
					}else if(bread){
						a=r.toText(true);
						if(paired){b=(r2==null ? new ByteBuilder(".") : r2.toText(true));}
					}else{
						throw new RuntimeException("Unsupported output format.");
					}
					
					a.append('\n');
					tsw1[mod].print(a);
					if(paired){
						b.append('\n');
						if(tsw2[i]!=null){tsw2[i].print(b);}
						else{tsw1[i].print(b);}
					}
					
					x++;
				}
			}
			
			cris.returnList(ln.id, readlist.isEmpty());
			
			//System.err.println("Waiting on a list...");
			ln=cris.nextList();
			readlist=ln.list;
		}
		
		//System.err.println("Returning a list... (final)");
		assert(readlist.isEmpty());
		cris.returnList(ln.id, readlist.isEmpty());
		
		
		for(TextStreamWriter tsw : tsw1){if(tsw!=null){tsw.poison();}}
		for(TextStreamWriter tsw : tsw2){if(tsw!=null){tsw.poison();}}
		ReadWrite.closeStream(cris);
		return x;
	}
	
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	public static int partitions=2;
	public static boolean fastq=false;
	public static boolean fasta=false;
	public static boolean bread=false;
	
}
