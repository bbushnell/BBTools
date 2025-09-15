package jgi;


import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import dna.Data;
import dna.Scaffold;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.LineParserS1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.LongList;
import tracker.ReadStats;

/**
 * 
 * Processes a sam file of mapped ESTs.
 * These ESTs may have been broken into smaller pieces for mapping,
 * and if so, are reassembled.
 * 
 * Produces a mapping statistics file.
 * 
 * @author Brian Bushnell
 * @date Sep 27, 2013
 *
 */
public class SamToEst {
	
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ReadWrite.USE_UNPIGZ=true;
		
		String est=null, stats=null, ref=null, sam=null;
		float fractionForAllCaptured=0.98f;
		
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
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("sam")){
				sam=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("stats")){
				stats=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("est")){
				est=b;
			}else if(a.equals("fraction")){
				fractionForAllCaptured=Float.parseFloat(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(sam==null && i==0 && Tools.looksLikeInputStream(arg)){
				sam=arg;
			}else if(stats==null && i==1 && !arg.contains("=")){
				stats=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		if(stats==null){stats="stdout";}
		SamToEst ste=new SamToEst(sam, stats, ref, est, fractionForAllCaptured);
		ste.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public SamToEst(String in_, String stats_, String ref_, String est_, float fractionForAll_){
		in=in_;
		stats=stats_;
		ref=ref_;
		estFile=est_;
		fractionForAll=fractionForAll_;
	}
	
	public void process(){
		HashMap<String, EST> table=new HashMap<String, EST>(initialSize);
		TextFile tf=new TextFile(in, true);
		String line=null;
		
		String program=null;
		String version=null;
		
		boolean bbmap=false;
		float bbversion=-1;
		LineParserS1 lp=new LineParserS1('\t');
		
		for(line=tf.nextLine(); line!=null && line.startsWith("@"); line=tf.nextLine()){
			lp.set(line);
			final String a=lp.parseString(1);
			
			if(a.equals("@SQ")){
				Scaffold sc=new Scaffold(lp);
//				assert(!table.containsKey(sc.name)) : "\nDuplicate scaffold name!\n"+sc+"\n\n"+table.get(sc.name);
//				table.put(sc.name, sc);
				refBases+=sc.length;
				refCount++;
			}else if(a.equals("@PG")){
				lp.set(line);
				for(int i=1, terms=lp.terms(); i<terms; i++) {
					String s=lp.parseString(i);
					if(s.startsWith("PN:")){
						String s2=s.substring(3);
						if(s2.equalsIgnoreCase("bbmap") || s2.startsWith("BBMap")){bbmap=true;}
						if(program==null){program=Data.forceIntern(s.substring(3));}
					}else if(s.startsWith("VN:")){
						if(bbmap && bbversion<0){bbversion=Float.parseFloat(s.substring(3));}
						if(version==null){version=Data.forceIntern(s.substring(3));}
					}
				}
			}else if(a.equals("@RG")){
				//Do nothing
			}else if(a.equals("@HD")){
				//Do nothing
			}else if(a.equals("@CO")){
				//Do nothing
			}else{
//				assert(false) : line;
			}
		}
		
		EST current=null;
		boolean err=false;
		for(; line!=null; line=tf.nextLine()){
			
			if(line.length()==0){
				
			}else if(line.charAt(0)=='@'){
				if(!err){
					outstream.println("Unexpected header line: "+line);
					outstream.println("This should not cause problems, and is probably due to concatenated sam files.\n" +
							"Supressing future unexpected header warnings.");
					err=true;
				}
				
				if(line.startsWith("@SQ")){
					String[] split=line.split("\t");
					Scaffold sc=new Scaffold(split);
//					if(!table.containsKey(sc.name)){
//						table.put(sc.name, sc);
//						refBases+=sc.length;
//						refCount++;
//					}
				}
			}else{
				
				SamLine sl=new SamLine(line);
				if(USE_SECONDARY || sl.primary()){
					
					if(sl.mapped() && sl.cigar!=null){
						String cigar=sl.cigar;
						if(cigar.contains("D") || cigar.contains("N")){
							int len=0;
							for(int i=0; i<cigar.length(); i++){
								char c=cigar.charAt(i);
								if(Tools.isDigit(c)){
									len=(len*10)+(c-'0');
								}else{
									if(c=='D' || c=='N'){
										introns.increment(len, 1);
									}
									len=0;
								}
							}
						}
					}
					
//					final Scaffold scaf=table.get(new String(sl.rname()));
//					assert(scaf!=null) : "Can't find "+new String(sl.rname());
//					final int a=Tools.max(sl.start(), 0);
//					final int b=Tools.min(sl.stop2(), scaf.length-1);
//					scaf.basehits+=(b-a+1);
					String name=sl.qname;
					int x=name.lastIndexOf('_');
					int part=1;
//					if(x>0){
					if(x>5 && name.charAt(x-5)=='_' && name.charAt(x-4)=='p' && name.charAt(x-3)=='a' && name.charAt(x-2)=='r' && name.charAt(x-1)=='t'){
						int partlen=name.length()-x-1;
						if(partlen>0 && partlen<6){
							int p2=0;
							for(int i=x+1; i<name.length(); i++){
								char c=name.charAt(i);
								int c2=c-'0';
								if(c2<0 || c2>9){
									p2=-1;
									break;
								}
							}
							if(p2>-1){
								part=p2;
								name=name.substring(0, x-5);
//								name=name.substring(0, x);
//								if(current!=null && !current.name.equals(name)){
//									//Special case test for sequences that already end with underscore number
//									if(name.length()>current.name.length()+1 && name.startsWith(current.name) && name.charAt(current.name.length())=='_'){
//										boolean specialCase=true;
//										for(int i=x+1; i<name.length(); i++){
//											char c=name.charAt(i);
//											int c2=c-'0';
//											if(c2<0 || c2>9){
//												specialCase=false;
//												break;
//											}
//										}
//										if(specialCase){name=current.name;}
//									}
//								}
							}else{
//								assert(false) : x+"\t"+p2+"\t"+name;
							}
						}else{
//							assert(false) : x+"\t"+name;
						}
					}else{
//						assert(false) : x+"\t"+name;
					}
					if(current==null || !current.name.equals(name)){
//						assert(part==1) : "Sam file must be in input order.  Run BBMap with the 'ordered' flag.\n"+part+"\n"+sl.qname;
						if(current!=null){addEst(current);}
						current=new EST(name);
					}
					current.add(sl);
				}
			}
		}
		if(current!=null){addEst(current);}
		tf.close();
		
		if(stats!=null){
			final TextStreamWriter tsw=new TextStreamWriter(stats, overwrite, false, false);
			tsw.start();
			
//			numRef:              786
//			numEst:            30985
//			EST-good:          30312 (   97.83%)
//			EST-best:          30312 (   97.83%)
//			EST-miss:            379 (    1.22%)
//			EST-zero:            294 (    0.95%)
			
//			tsw.println("EST-good:\t"+good+"\t"++"");
//			tsw.println("EST-best:\t"+best+"\t"++"");
//			tsw.println("EST-miss:\t"+miss+"\t"++"");
//			tsw.println("EST-zero:\t"+zero+"\t"++"");
			
			boolean oldStyle=false;

			if(oldStyle){
				tsw.println("ref:\t"+ref);
				tsw.println("est:\t"+estFile);
				tsw.println("sam:\t"+in);

				tsw.println("numRef:\t"+refCount+"\t"+refBases);
				tsw.println("numEst:\t"+estCount+"\t"+estBases);
				tsw.println("type\t#ests\t%ests\t#bases\t%bases");
			}else{

				tsw.println("ref_file="+ref);
				tsw.println("est_file="+estFile);
				tsw.println("sam_file="+in);

				tsw.println("n_ref_scaffolds="+refCount);
				tsw.println("n_ref_bases="+refBases);
				tsw.println("n_est="+estCount);
				tsw.println("n_est_bases="+estBases);
				tsw.println("type\tn_est\tpct_est\tn_bases\tpct_bases");
			}
			
			double multE=100.0/estCount;
			double multB=100.0/estBases;

			double allBasesPct=multE*allBasesMapped;
			double mostBasesPct=multE*mostBasesMapped;
			double someBasesPct=multE*someBasesMapped;
			double noBasesPct=multE*noBasesMapped;
			double multiScaffoldPct=multE*multiScaffold;
			
			double allBasesPctB=multB*allBasesMappedB;
			double mostBasesPctB=multB*mostBasesMappedB;
			double someBasesPctB=multB*someBasesMappedB;
			double noBasesPctB=multB*noBasesMappedB;
			double multiScaffoldPctB=multB*multiScaffoldB;
			
			int min=0, max=0, median=0;
			long sum=0, count=0;
			for(int i=minIntron; i<introns.size; i++){
				long x=introns.get(i);
				if(x>0){
					if(min==0){min=i;}
					max=i;
					sum+=(i*x);
					count+=x;
				}
			}
			if(count>0){ //If there are any introns
				long half=(count+1)/2; //50th percentile of number of introns
				assert(half<=count);
				long count2=0; //Current sum of length
				for(int i=0; count2<half; i++){
					long x=introns.get(i);
					if(x>0){
						count2+=x;
						median=i;
					}
				}
			}
			
			tsw.println("all:\t"+allBasesMapped+"\t"+Tools.format("%.4f%%",allBasesPct)+"\t"+allBasesMappedB+"\t"+Tools.format("%.4f%%",allBasesPctB));
			tsw.println("most:\t"+mostBasesMapped+"\t"+Tools.format("%.4f%%",mostBasesPct)+"\t"+mostBasesMappedB+"\t"+Tools.format("%.4f%%",mostBasesPctB));
			tsw.println("some:\t"+someBasesMapped+"\t"+Tools.format("%.4f%%",someBasesPct)+"\t"+someBasesMappedB+"\t"+Tools.format("%.4f%%",someBasesPctB));
			tsw.println("zero:\t"+noBasesMapped+"\t"+Tools.format("%.4f%%",noBasesPct)+"\t"+noBasesMappedB+"\t"+Tools.format("%.4f%%",noBasesPctB));
			tsw.println("multi:\t"+multiScaffold+"\t"+Tools.format("%.4f%%",multiScaffoldPct)+"\t"+multiScaffoldB+"\t"+Tools.format("%.4f%%",multiScaffoldPctB));
//			tsw.println("numIntrons:\t"+count);
//			tsw.println("minIntron:\t"+min);
//			tsw.println("maxIntron:\t"+max);
//			tsw.println("medIntron:\t"+median);
//			tsw.println("avgIntron:\t"+(long)(sum/(double)(Tools.max(count,1))));
			tsw.println("introns\tmin\tmax\tmedian\taverage");
			tsw.println(count+"\t"+min+"\t"+max+"\t"+median+"\t"+Tools.format("%.1f", (sum/(double)(Tools.max(count,1)))));

			tsw.poisonAndWait();
		}
	}

	private void addEst(EST est){
//		outstream.println("\n"+est);
		estCount++;
		partCount+=est.parts;
		estBases+=est.length;
		estBasesMapped+=est.mappedLength;
		partCountMapped+=est.mappedParts;
		
		for(int i=0; i<est.msdicn.length; i++){
			msdicnOverall[i]+=est.msdicn[i];
		}
		
		if(est.scafnames.size()>1){
			multiScaffold++;
			multiScaffoldB+=est.length;
		}
		
		if(est.mappedParts==est.parts){
//			outstream.print("A");
			allPartsMapped++;
		}else if(est.mappedParts>=Tools.max(1, est.parts/2)){
//			outstream.print("B");
			mostPartsMapped++;
		}else if(est.mappedParts>0){
//			outstream.print("C");
			somePartsMapped++;
		}else{
//			outstream.print("D");
			noPartsMapped++;
		}
		
		int match=est.match();
		if(match>=(est.length*fractionForAll)){
//			outstream.print("E");
			allBasesMapped++;
			allBasesMappedB+=est.length;
		}else if(match>=est.length/2){
//			outstream.print("F");
			mostBasesMapped++;
			mostBasesMappedB+=est.length;
		}else if(match>0){
//			outstream.print("G");
			someBasesMapped++;
			someBasesMappedB+=est.length;
		}else{
//			outstream.print("H");
			noBasesMapped++;
			noBasesMappedB+=est.length;
		}
	}
	
	public final float fractionForAll;
	public final String in, stats, ref, estFile;
	
	public long refBases=0;
	public long estBases=0;
	public long estBasesMapped=0;

	public long refCount=0;
	public long estCount=0;
	public long partCount=0;
	public long partCountMapped=0;
	
	public long good=0, best=0, miss=0, zero=0;
	public long multiScaffold=0, multiScaffoldB=0;
	public long allPartsMapped=0, mostPartsMapped=0, somePartsMapped=0, noPartsMapped=0;
	public long allBasesMapped=0, mostBasesMapped=0, someBasesMapped=0, noBasesMapped=0;
	public long allBasesMappedB=0, mostBasesMappedB=0, someBasesMappedB=0, noBasesMappedB=0;
	public long[] msdicnOverall=new long[6];
	public LongList introns=new LongList(1);
	
	public int initialSize=4096;
	public boolean ADD_FROM_REF=true;
	public boolean USE_SECONDARY=false;
	public static int minIntron=10;
	public static boolean overwrite=true;
	public static boolean append=false;
//	public HashMap<String, EST> //Only needed if sam file is unordered.
	
	static PrintStream outstream=System.err;
	
	public static class EST{
		
		public EST(String name_){
			name=name_;
			outstream.println("New EST: "+name);
		}
		
		public void add(SamLine sl){
			outstream.println("Adding samline "+sl.qname+" to EST "+name);
			parts++;
//			length+=sl.seq.length();
			length+=sl.seq.length;
			if(sl.mapped()){
//				mappedLength+=sl.seq.length();
				mappedLength+=sl.seq.length;
				mappedParts++;
				if(sl.cigar!=null){
					String matchTag=sl.matchTag();
					
					int[] temp;
					if(matchTag==null){
						temp=SamLine.cigarToMsdic(sl.cigar);
					}else{
						temp=Read.matchToMsdicn(matchTag.getBytes());
					}
					for(int i=0; i<temp.length; i++){
						msdicn[i]+=temp[i];
					}
				}
				if(sl.rname()!=null){
					scafnames.add(new String(sl.rname()));
				}
			}
		}
		
		public int match(){return msdicn[0];}
		
		@Override
		public String toString(){
			StringBuilder sb=new StringBuilder();
			sb.append(name).append('\t');
			sb.append(length).append('\t');
			sb.append(mappedLength).append('\t');
			sb.append(parts).append('\t');
			sb.append(mappedParts).append('\t');
			sb.append(Arrays.toString(msdicn)).append('\t');
			sb.append(scafnames).append('\t');
			return sb.toString();
		}
		
		final String name;
		int length=0, mappedLength=0;
		int parts=0, mappedParts=0;
		HashSet<String> scafnames=new HashSet<String>(4);
		
		int[] msdicn=KillSwitch.allocInt1D(6);
		
	}
	
}
