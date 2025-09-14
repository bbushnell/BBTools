package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.LongHashSet;
import tax.TaxFilter;

/**
 * From a 3+ column text file of {node 1, node 2, dist, optionally sizeratio and others},
 * makes a minimal representative set by retaining nodes
 * such that all original nodes are within a minimum distance
 * of at least one representative node.
 * Singleton nodes will only be included if they are
 * represented by a self-edge.
 * 
 * @author Brian Bushnell
 * @date October 26, 2017
 *
 */
public class RepresentativeSet {
	
	public static void main(String[] args){
		Timer t=new Timer();
		RepresentativeSet x=new RepresentativeSet(args);
		x.process(t);
		x.printResults(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public RepresentativeSet(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		boolean taxFlag=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("thresh") || a.equals("threshold") || a.equals("minid") || a.equals("minani") || a.equals("id") || a.equals("ani")){
				threshold=Float.parseFloat(b);
			}else if(a.equals("ratio") || a.equals("minratio") || a.equals("sizeratio") || a.equals("sr")){
				minRatio=Float.parseFloat(b);
			}
			
			else if(a.equals("maxsize")){
				maxSize=Parse.parseKMG(b);
			}else if(a.equals("minsize")){
				minSize=Parse.parseKMG(b);
			}
			else if(a.equals("maxbases") || a.equals("maxbp")){
				maxBases=Parse.parseKMG(b);
			}else if(a.equals("minbases") || a.equals("minbp")){
				minBases=Parse.parseKMG(b);
			}
			
			else if(a.equals("invertratio") || a.equals("ir")){
				invertRatio=Parse.parseBoolean(b);
			}else if(a.equals("printsize")){
				printSize=Parse.parseBoolean(b);
			}else if(a.equals("printclusters") || a.equals("cluster")){
				printClusters=Parse.parseBoolean(b);
			}else if(a.equals("printheader") || a.equals("header")){
				printHeader=Parse.parseBoolean(b);
			}
			
//			else if(a.equals("taxfilter")){
//				foo=Parse.parseBoolean(b);
//			}
			
			else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(TaxFilter.validArgument(a)){
				taxFlag=true;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		
		if(taxFlag){
			filter=TaxFilter.makeFilter(args);
			if(filter.isEmpty()){filter=null;}
		}
	}
	
	void process(Timer t){
		
		HashMap<Long, Node> map=load();
		LongHashSet set=new LongHashSet();

		ByteStreamWriter bsw=null;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}

		ArrayList<Node> list=new ArrayList<Node>(map.size());
		ArrayList<Node> singletons=new ArrayList<Node>(map.size());
		{//Isolate singletons to accelerate sorting.  This is not really necessary.
			ArrayList<Node> list0=new ArrayList<Node>(map.size());
			list0.addAll(map.values());
			for(Node n : list0){
				sizeProcessed+=n.size;
				basesProcessed+=n.bases;
				if(n.edges==null || n.edges.isEmpty()){
					singletons.add(n);
					nodesValid++;
					sizeValid+=n.size;
					basesValid+=n.bases;
					n.used=true;
				}else{
					list.add(n);
				}
			}
		}
		
		Collections.sort(list);
		Collections.reverse(list);
//		assert(list.size()==0 || list.get(0).sum>=list.get(list.size()-1).sum) : "\n"+list+"\n"; //This is incorrect because the comparator changed
		
		for(Node n : list){
			assert(!n.used);
			boolean ok=true;
			if(n.edges!=null){
				for(Edge e : n.edges){
					if(set.contains(e.b)){
						ok=false;
						break;
					}
				}
			}
			if(ok){
				set.add(n.id);
				n.used=true;
			}
		}
		
		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}
		
		if(printHeader){
			StringBuilder header=new StringBuilder();
			header.append("#Representative");
			if(printSize){header.append("\tSize");}
			if(printClusters){header.append("\tNodeCount\tNodes");}
			if(bsw!=null){bsw.println(header);}
			if(bswInvalid!=null){bswInvalid.println(header);}
		}
		
		for(Node n : list){
			if(n.used){
				nodesValid++;
				sizeValid+=n.size;
				basesValid+=n.bases;
				if(bsw!=null){
					bsw.print(n.id);
					if(printSize){
						bsw.print('\t').print(n.size);
					}
					if(printClusters){
						bsw.print('\t').print(n.edges.size());
						bsw.print('\t').print(n.edges.get(0).b);
						for(int i=1; i<n.edges.size(); i++){
							bsw.print(',').print(n.edges.get(i).b);
						}
					}
					bsw.println();
				}
			}else{
				if(bswInvalid!=null){
					bswInvalid.print(n.id);
					if(printSize){
						bswInvalid.print('\t').print(n.size);
					}
					if(printClusters){
						bswInvalid.print('\t').print(n.edges.size());
						bswInvalid.print('\t').print(n.edges.get(0).b);
						for(int i=1; i<n.edges.size(); i++){
							bswInvalid.print(',').print(n.edges.get(i).b);
						}
					}
					bswInvalid.println();
				}
			}
		}
		for(Node n : singletons){
			if(bsw!=null){
				bsw.print(n.id);
				if(printSize){
					bsw.print('\t').print(n.size);
				}
				if(printClusters && n.edges!=null){
					bsw.print('\t').print(n.edges.size());
					bsw.print('\t').print(n.edges.get(0).b);
					for(int i=1; i<n.edges.size(); i++){
						bsw.print(',').print(n.edges.get(i).b);
					}
				}
				bsw.println();
			}
		}
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
	}
	
	private HashMap<Long, Node> load(){
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		
		HashMap<Long, Node> map=new HashMap<Long, Node>();
		
		byte[] line=bf.nextLine();
		if(line!=null && (Tools.startsWith(line, "Query\t") || line[0]=='#')){
			line=bf.nextLine();
		}
		
		LongHashSet ignored=new LongHashSet();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				Edge e=new Edge(line);
				boolean pass=false;
				if(e.sizeA<1 || (e.sizeA>=minSize && (maxSize<1 || e.sizeA<=maxSize))){//Size filter
					if(e.basesA<1 || (e.basesA>=minBases && (maxBases<1 || e.basesA<=maxBases))){//Bases filter
						pass=true;
						if(filter==null || (filter.passesFilter((int)e.a) && filter.passesFilter((int)e.b))){//Tax filter
							if(e.dist>=threshold && e.ratio()>=minRatio){
								Node n=map.get(e.a);
								if(n==null){
									nodesProcessed++;
									n=new Node(e.a, e.sizeA, e.basesA);
									map.put(e.a, n);
								}
								if(e.a!=e.b){n.add(e);}
							}
						}
					}
				}
				if(!pass && !ignored.contains(e.a)){
					ignored.add(e.a);
					nodesIgnored++;
					sizeIgnored+=e.sizeA;
					basesIgnored+=e.basesA;
				}
			}
			line=bf.nextLine();
		}
		
		errorState|=bf.close();
		return map;
	}
	
	private void printResults(Timer t){
		t.stop();

		long np2=nodesProcessed+nodesIgnored;
		double rpnano=linesProcessed/(double)(t.elapsed);
		double npnano=np2/(double)(t.elapsed);
		double bpnano=bytesProcessed/(double)(t.elapsed);

		String rpstring=(linesProcessed<100000 ? ""+linesProcessed : linesProcessed<100000000 ? (linesProcessed/1000)+"k" : (linesProcessed/1000000)+"m");
		String npstring=(np2<100000 ? ""+np2 : np2<100000000 ? (np2/1000)+"k" : (np2/1000000)+"m");
		String bpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(npstring.length()<8){npstring=" "+npstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Lines Processed:    "+rpstring+" \t"+Tools.format("%.2fk lines/sec", rpnano*1000000));
		outstream.println("Nodes Processed:    "+npstring+" \t"+Tools.format("%.2fk lines/sec", npnano*1000000));
		outstream.println("Bytes Processed:    "+bpstring+" \t"+Tools.format("%.2fm bytes/sec", bpnano*1000));
		
		outstream.println();
		outstream.println("Valid Nodes:       \t"+nodesValid);
		outstream.println("Invalid Nodes:     \t"+(nodesProcessed-nodesValid));
		outstream.println("Ignored Nodes:     \t"+nodesIgnored);
		outstream.println();
		outstream.println("Valid Size:        \t"+sizeValid);
		outstream.println("Invalid Size:      \t"+(sizeProcessed-sizeValid));
		outstream.println("Ignored Size:      \t"+sizeIgnored);
		outstream.println();
		outstream.println("Valid Bases:        \t"+basesValid);
		outstream.println("Invalid Bases:      \t"+(basesProcessed-basesValid));
		outstream.println("Ignored Bases:      \t"+basesIgnored);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	
	private class Node implements Comparable<Node> {
		
		Node(long id_, long size_, long bases_){
			id=id_;
			size=Tools.max(2, size_);
			bases=Tools.max(2, bases_);
//			assert(bases>2) : bases+", "+size;
		}
		
		void add(Edge e){
			if(edges==null){edges=new ArrayList<Edge>();}
			edges.add(e);
			sum+=e.dist-threshold+0.000001;
		}
		
		@Override
		public int compareTo(Node b) {
			double ca=sum+0.25*sum*Math.log(size);
			double cb=b.sum+0.25*b.sum*Math.log(b.size);
			if(ca>cb){return 1;}
			if(ca==cb){return 0;}
			return -1;
//			if(sum>b.sum){return 1;}
//			if(sum==b.sum){return 0;}
//			return -1;
		}
		
		@Override
		public String toString(){
			StringBuilder sb=new StringBuilder();
			sb.append("Node ").append(id);
			sb.append(", size ").append(size);
			sb.append(", sum ").append(sum);
			sb.append(", edges: {");
			char c=0;
			for(Edge e : edges){
				if(c>0){sb.append(c);sb.append(' ');}
				else{c=',';}
				sb.append(e);
			}
			sb.append('}').append('\n');
			return sb.toString();
		}
		
		final long id;
		final long size;
		final long bases;
		ArrayList<Edge> edges;
		boolean used=false;
		double sum;
		
	}
	
	private class Edge {
		
		Edge(byte[] line){
			String[] split=new String(line).split("\t+");
			try {
				a=Long.parseLong(split[0]);
				b=Long.parseLong(split[1]);
				dist=Float.parseFloat(split[2]);
			} catch (NumberFormatException e) {
				System.err.println(new String(line)+"\n"+Arrays.toString(split));
				e.printStackTrace();
				KillSwitch.kill();
				throw new RuntimeException();
			}
			long sizeA_=1, sizeB_=1;
			if(split.length>4){
				try {
					sizeA_=Long.parseLong(split[3]);
					sizeB_=Long.parseLong(split[4]);
				} catch (NumberFormatException e) {}
			}
			sizeA=sizeA_;
			sizeB=sizeB_;
			long basesA_=1, basesB_=1;
			if(split.length>6){
				try {
					basesA_=Long.parseLong(split[5]);
					basesB_=Long.parseLong(split[6]);
				} catch (NumberFormatException e) {}
			}
			basesA=basesA_;
			basesB=basesB_;
//			assert(basesA>2) : new String(line);
		}
		
		Edge(long a_, long b_, double dist_, long sizeA_, long sizeB_, long basesA_, long basesB_){
			a=a_;
			b=b_;
			dist=dist_;
			sizeA=sizeA_;
			sizeB=sizeB_;
			basesA=basesA_;
			basesB=basesB_;
		}
		
		double ratio(){
			double ratio=Tools.max(1, sizeA)/(float)Tools.max(1, sizeB);
			return (invertRatio && ratio>1 ? 1/ratio : ratio);
		}
		
		@Override
		public String toString(){
			ByteBuilder sb=new ByteBuilder();
			sb.append('(').append(a).append(',').append(b).append(',').append(dist, 3).append(')');
			return sb.toString();
		}
		
		final long a;
		final long b;
		final double dist;
		final long sizeA;
		final long sizeB;
		final long basesA;
		final long basesB;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/

	private long nodesProcessed=0;
	private long nodesValid=0;
	private long sizeValid=0;
	private long basesValid=0;
	
	private long linesProcessed=0;
	private long bytesProcessed=0;
	private long sizeProcessed=0;
	private long basesProcessed=0;
	
	private long nodesIgnored=0;
	private long sizeIgnored=0;
	private long basesIgnored=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	private double threshold=0;
	private double minRatio=0;
	
	/** Ignore nodes over this estimated size in unique kmers */
	private long maxSize=-1;
	/** Ignore nodes under this estimated size in unique kmers */
	private long minSize=0;

	/** Ignore nodes over this estimated size in total bp */
	private long maxBases=-1;
	/** Ignore nodes under this estimated size in total bp */
	private long minBases=0;
	
	boolean invertRatio=false;
	boolean printSize=true;
	boolean printHeader=true;
	boolean printClusters=true;
	
	private TaxFilter filter;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private final FileFormat ffoutInvalid;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
