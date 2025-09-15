package tax;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.IntList;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class FindAncestor {
	
	public static void main(String[] args){
		Timer t=new Timer();
		FindAncestor x=new FindAncestor(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public FindAncestor(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
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
		
		if(giTableFile!=null){
			outstream.println("Loading gi table.");
			GiToTaxid.initialize(giTableFile);
		}
		if(taxTreeFile!=null){
			tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, true);
			assert(tree.nameMap!=null);
		}else{
			tree=null;
			throw new RuntimeException("No tree specified.");
		}
		lifeNode=tree.getNodeByName("life");
	}
	
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		
		bsw.print("#Name\tAncestor\tMajority\tTaxonomy...\n".getBytes());
		
		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}
		
//		final HashArray1D counts=countTable ? new HashArray1D(256000, true) : null;
		final IntList giList=new IntList();
		final IntList tidList=new IntList();
		final IntList traversal=new IntList();
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				
				giList.clear();
				tidList.clear();
				traversal.clear();
				
				final int giCount=getGiNumbers(line, giList, ',');
				final int ncbiCount=getTaxidNumbers(giList, tidList);
				
				taxaCounted+=giCount;
				taxaValid+=ncbiCount;
				final boolean valid=(ncbiCount>0);
				
				if(valid){
					linesValid++;
					int ancestor=findAncestor(tidList);
					int majority=findMajority(tidList);
					
					for(int i=0; i<line.length && line[i]!='\t'; i++){
						bb.append(line[i]);
					}
					bb.tab();
					bb.append(ancestor);
					bb.tab();
					bb.append(majority);
					bb.tab();
					
					fillTraversal(majority, traversal, true);
					writeTraversal(traversal, bb);
					bb.nl();
					
					for(int i=0; i<tidList.size; i++){
						final int id=tidList.get(i);
						fillTraversal(id, traversal, true);
						writeTraversal(traversal, bb);
						bb.nl();
					}
					
					bsw.print(bb.toBytes());
					bb.clear();
				}else{
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
			}
			line=bf.nextLine();
		}
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Valid Lines:       \t"+linesValid);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesValid));
//		if(counts!=null){
//			outstream.println("Unique Taxa:       \t"+taxaCounted);
//		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void fillTraversal(int id, IntList traversal, boolean addLife){
		traversal.clear();
		for(TaxNode node=tree.getNode(id); node!=null && node!=lifeNode; node=tree.getNode(node.pid)){
			traversal.add(node.id);
		}
		if(addLife || traversal.size==0){traversal.add(lifeNode.id);}
	}
	
	private void writeTraversal(IntList traversal, ByteBuilder bb){
		for(int i=traversal.size-1; i>=0; i--){
			final int id=traversal.get(i);
			if(id>=0){
				TaxNode tn=tree.getNode(id);
				//			bb.append(tn.level+"_"+tn.name);
				bb.append(/*tn.level+"_"+*/tn.name);
				if(i>0){bb.tab();}
			}
		}
	}
	
	private int getGiNumbers(final byte[] line, final IntList list, final char delimiter){
		int i=0;
		
		//Skip name
		while(i<line.length && line[i]!='\t'){i++;}
		
		//Skip whitespaces
		while(i<line.length && Character.isWhitespace(line[i])){i++;}
		
		while(i<line.length){
			while(i<line.length && line[i]==delimiter){i++;}
			int start=i;
			while(i<line.length && line[i]!=delimiter){i++;}
			final int stop=i;
			if(Tools.startsWith(line, prefix, start)){start+=3;}
			assert(start<stop) : "Badly formatted line at "+start+":\n"+new String(line);
//			System.err.println(start+","+stop+",'"+new String(line).substring(start, stop)+"'");
			if(start<stop){
				final int number=Parse.parseInt(line, start, stop);
				list.add(number);
			}
		}
		return list.size;
	}
	
	private static int getTaxidNumbers(final IntList giList, final IntList ncbiList){
		final int size=giList.size;
		for(int i=0; i<size; i++){
			final int gi=giList.get(i);
			final int ncbi=GiToTaxid.getID(gi);
//			System.err.println(gi+" -> "+ncbi);
			if(ncbi>=0){ncbiList.add(ncbi);}
		}
		return ncbiList.size;
	}
	
	private int findAncestor(IntList list){
		return findAncestor(tree, list);
	}
	
	public static int findAncestor(TaxTree tree, IntList list){
		if(list.size<1){
			assert(false);
			return -1;
		}
		int ancestor=list.get(0);
		for(int i=1; i<list.size && ancestor>-1; i++){
			final int id=list.get(i);
//			System.err.println(ancestor+"+"+id+" -> "+tree.commonAncestor(ancestor, id));
			int x=tree.commonAncestor(ancestor, id);
			if(x>-1){
				ancestor=x;
			}
		}
//		System.err.println("Ancestor node: "+tree.getNode(ancestor));
//		System.err.println(list+" -> "+ancestor);
//		if(ancestor<0){ancestor=lifeNode.id;}
		return ancestor;
	}
	
	private int findMajority(IntList list){
		if(list.size<3){return findAncestor(list);}
		final int majority=list.size/2+1;
//		System.err.println("Majority: "+majority);
		
		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
//			System.err.println("Found node "+tn);
			assert(tn!=null) : "No node for id "+id;
			if(tn!=null){
				tree.percolateUp(tn, 1);
			}
		}
		
		TaxNode best=lifeNode;
		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
			while(tn!=null && tn!=lifeNode){
				if(tn.countSum>=majority && tn.levelExtended<best.levelExtended){
					best=tn;
					break;
				}
				tn=tree.getNode(tn.pid);
			}
		}
		
//		System.err.println("Best node: "+best);
		
		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
			if(tn!=null){
				tree.percolateUp(tn, -1);
			}
		}
		
		return best.id;
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String outInvalid=null;

	private String giTableFile=null;
	private String taxTreeFile=null;
	
	private final TaxTree tree;
	
	private final TaxNode lifeNode;
	
	/*--------------------------------------------------------------*/
	
	private long taxaCounted=0;
	private long taxaValid=0;
	private long linesProcessed=0;
	private long linesValid=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;

//	private boolean prefix=false;
	private boolean countTable=true;
//	private boolean keepInvalidSequence=false;
	
	private final String prefix="gi|";
	
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
