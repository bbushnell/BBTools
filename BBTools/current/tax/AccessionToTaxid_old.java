package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicLongArray;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.HashBuffer;
import kmer.KmerTableSet;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.StringNum;

/**
 * @author Brian Bushnell
 * @date December 16, 2016
 *
 */
public class AccessionToTaxid_old {
	
	public static void load(String files){
		final boolean oldBf2=ByteFile.FORCE_MODE_BF2;
		final boolean oldBf1=ByteFile.FORCE_MODE_BF1;
		final boolean oldUnpigz=ReadWrite.USE_UNPIGZ;
		final boolean oldGunzip=ReadWrite.USE_UNPIGZ;
		
		main(new String[] {"in="+files, "unpigz="+ReadWrite.USE_UNPIGZ, "gunzip="+ReadWrite.USE_GUNZIP});

		ByteFile.FORCE_MODE_BF2=oldBf2;
		ByteFile.FORCE_MODE_BF1=oldBf1;
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		ReadWrite.USE_UNPIGZ=oldGunzip;
	}
	
	public static void main(String[] args){
		Timer t=new Timer();
		AccessionToTaxid_old x=new AccessionToTaxid_old(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public AccessionToTaxid_old(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_UNPIGZ=true;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("stripunderscore")){
//				STRIP_UNDERSCORE=Parse.parseBoolean(b);
				assert(false) : "stripunderscore is disabled.";
			}else if(a.equals("usetables")){
//				USE_TABLES=Parse.parseBoolean(b);
			}else if(a.equals("usetables")){
//				USE_TABLES=Parse.parseBoolean(b);
			}else if(a.equals("skipparse")){
				skipParse=Parse.parseBoolean(b);
			}else if(a.equals("skiphash")){
				skipHash=Parse.parseBoolean(b);
			}else if(a.equals("prealloc")){ 
				if(b==null || Character.isLetter(b.charAt(0))){
					if(Parse.parseBoolean(b)){
						prealloc=0.78f;
					}else{
						prealloc=0;
					}
				}else{
					prealloc=Float.parseFloat(b);
				}
			}else if(a.equals("maxpigzprocesses")){
				maxPigzProcesses=Integer.parseInt(b);
			}else if(a.equals("in")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] temp=b.split(",");
				for(String s : temp){in.add(s);}
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(b==null){
				if(new File(arg).exists()){
					in.add(arg);
				}
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;

//			out=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null || in.size()==0){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(ReadWrite.USE_UNPIGZ && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

//		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
//		if(!Tools.testOutputFiles(overwrite, false, false, out)){
//			outstream.println((out==null)+", "+out);
//			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
//		}

		{//Reorder by size, ascending
			ArrayList<StringNum> list=new ArrayList<StringNum>();
			for(String s : in){
				list.add(new StringNum(s, new File(s).length()));
			}
			Collections.sort(list);
			in.clear();
			for(StringNum sn : list){
				in.add(sn.s);
			}
		}
		
//		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, false, false);
		ffin=new FileFormat[in.size()];
		
		/* Note */
		/* Java 1.7 works fine here (54 seconds skipping parsing). */
		/* Java 1.8 has immense speed-downs if pigz is used (80-100s normally, >1000s with unpigz). */
		/* Java 1.8_144 is unpredictable and incredibly slow (80-900s normally, 500-1800 with unpigz) */
		
		int processes=0;
		for(int i=0; i<in.size(); i++){
			String s=in.get(i);
			if(!new File(s).exists()){
				if(s.startsWith("shrunk.") && new File(s.substring(7)).exists()){
					s=s.substring(7);
				}
			}
			FileFormat ff=FileFormat.testInput(s, FileFormat.TXT, null, true, false);
			if(ff.gzip() && processes>maxPigzProcesses){
				processes++;
//				if(processes>maxPigzProcesses){
					ff=FileFormat.testInput(s, FileFormat.TXT, null, false, false);
//				}
			}
			ffin[i]=ff;
		}
	}
	
	@SuppressWarnings("unchecked")
	void process(Timer t){

//		if(USE_MAPS){
			assert(maps==null);
			maps=new HashMap[128];
			for(int i=0; i<maps.length; i++){
				maps[i]=new HashMap<String, Integer>();
			}
//		}

		assert(tables==null);
		if(USE_TABLES){
			tables=new KmerTableSet(new String[] {"ways=31",("prealloc="+(prealloc>0 ? prealloc : "f"))}, 12);
			tables.allocateTables();
		}
		
		if(ffin.length>4){//Addresses a multithreaded read bug in Java
			FileFormat[] ffa1=Arrays.copyOf(ffin, 2);
			FileFormat[] ffa2=Arrays.copyOfRange(ffin, 2, ffin.length);
			spawnThreads(ffa1);
			spawnThreads(ffa2);
		}else{
			spawnThreads(ffin);
		}
		
		//Do anything necessary after processing
		System.gc();
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Valid Lines:       \t"+linesValid);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesValid));

		if(lengthCounts!=null){
			outstream.println();
			outstream.println("Length counts:");

			for(int i=0; i<lengthCounts.length(); i++){
				long count=lengthCounts.get(i);
				if(count>0){outstream.println(i+"\t"+count);}
			}
		}

		if(symbolCounts!=null){
			outstream.println();
			outstream.println("Symbols:");
			
			String comma="";
			for(int i=0; i<symbolCounts.length(); i++){
				long count=symbolCounts.get(i);
				if(count>0){
					outstream.print(comma+i);
					comma=",";
				}
			}
		}

		if(counts_underscore!=null){
			outstream.println();
			outstream.println("Length_underscore counts:");

			for(int i=0; i<counts_underscore.length(); i++){
				long count=counts_underscore.get(i);
				if(count>0){outstream.println(i+"\t"+count);}
			}
		}

		if(counts_underscore2!=null){
			outstream.println();
			outstream.println("Length_underscore2 counts:");

			for(int i=0; i<counts_underscore2.length(); i++){
				long count=counts_underscore2.get(i);
				if(count>0){outstream.println(i+"\t"+count);}
			}
		}
		outstream.println();
		Shared.printMemory();
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		
		LOADED=true;
	}
	
	/** Spawn process threads */
	private void spawnThreads(FileFormat[] ffa){
		
		//Do anything necessary prior to processing
		
		//Fill a list with ProcessThreads
		ArrayList<HashThread> alht=new ArrayList<HashThread>(ffa.length);
		for(FileFormat ff : ffa){
			if(ff!=null){
				System.err.println("Loading "+ff.name());
				alht.add(new HashThread(ff));
			}
		}
		
		//Start the threads
		for(HashThread pt : alht){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		for(HashThread pt : alht){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			linesProcessed+=pt.linesProcessedT;
			linesValid+=pt.linesValidT;
			bytesProcessed+=pt.bytesProcessedT;

			accumulate(lengthCounts, pt.lengthCountsT);
			accumulate(symbolCounts, pt.symbolCountsT);
			accumulate(counts_underscore, pt.counts_underscoreT);
			accumulate(counts_underscore2, pt.counts_underscore2T);
			
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	private static void accumulate(AtomicLongArray a, long[] b){
		if(a==null || b==null){return;}
		for(int i=0; i<b.length; i++){
			a.getAndAdd(i, b[i]);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	public static int get(String accession){
		if(accession==null){return -1;}
//		if(STRIP_UNDERSCORE){
//			accession=accession.replaceAll("[_-]", "");
//		}
		
		int dot=accession.indexOf('.');
		int len=(dot<0 ? accession.length() : dot);
		
		
		
		if(USE_TABLES){
			if(AnalyzeAccession.codeMap!=null){
//				if(dot>AnalyzeAccession.longestPattern){return false;}
				final long number=AnalyzeAccession.digitize(accession);
				if(number>=0){
					int value=tables.getCount(number);
					return value<0 ? -1 : value;
				}
			}else if(len<=12){
				long number=hash(accession);

				int value=tables.getCount(number);
				return value<1 ? -1 : value;
			}
		}
		
		if(dot>-1){accession=accession.substring(0,  dot);}
		if(accession.length()<1){return -1;}
		int way=accession.charAt(0);
		Integer value=maps[way].get(accession);
		return value==null ? -1 : value.intValue();
	}
	
	public static boolean isValidAccession(String s){
		if(s==null || s.length()<4){return false;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if((c>='0' && c<='9') || (c>='A' && c<='Z') /*|| (c>='a' && c<='z')*/ || c=='.' || c=='_' || c=='-'){
				//do nothing
			}else{
				return false;
			}
		}
		return true;
	}
	
	static long hash(String accession){
		long number=0;
		for(int i=0, max=accession.length(); i<max; i++){
			long c=accession.charAt(i);
			if(c=='.'){break;}
			if(c>='0' && c<='9'){c=c-'0';}
			else if(c>='A' && c<='Z'){c=c+offset;}
			else if(c=='_' || c=='-'){c=10;}//Collision, but should be OK
			else if(c>='a' && c<='z'){c=c+offsetLower;}
			else{
				assert(false) : accession;
			}
			number=(number*37)+c;
		}
		return number;
	}
	
	static long hash(final byte[] line, final int limit){
		long number=0;
		for(int i=0; i<limit; i++){
			long c=line[i];
			if(c=='.'){break;}
			if(c>='0' && c<='9'){c=c-'0';}
			else if(c>='A' && c<='Z'){c=c+offset;}
			else if(c=='_' || c=='-'){c=10;}//Collision, but should be OK
			else if(c>='a' && c<='z'){c=c+offsetLower;}
			else{
				assert(false) : new String(line);
			}
			number=(number*37)+c;
		}
		return number;
	}
	
	public static int parseLineToTaxid(final byte[] line, final byte delimiter){
		int a=0, b=0;
		
		final int ncbi;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
//		assert(b>a) : "Missing field 1: "+new String(line);
		assert(b>=a) : "Missing field 1: "+new String(line)+"\n"+a+", "+b;
		//accession2=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		ncbi=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		return ncbi;
	}
	
	/*--------------------------------------------------------------*/
	
	public static class HashThread extends Thread {
		
		@SuppressWarnings("unchecked")
		public HashThread(FileFormat ff_){
//			if(USE_MAPS){
				mapsT=new HashMap[128];
				for(int i=0; i<mapsT.length; i++){
					mapsT[i]=new HashMap<String, Integer>();
				}
//			}
			if(USE_TABLES){
				table=new HashBuffer(tables.tables(), 1000, 31, true, true);
			}
			ff=ff_;
		}
		
		@Override
		public void run(){
			
			ByteFile bf=ByteFile.makeByteFile(ff, 1);
			
			byte[] line=bf.nextLine();
			while(line!=null && Tools.startsWith(line, "accession")){line=bf.nextLine();}
			
			while(line!=null){
				if(line.length>0){
					linesProcessedT++;
					bytesProcessedT+=line.length;

//					final boolean valid=(!Tools.startsWith(line, "accession\t")) & !skipParse;
					final boolean valid=(!Tools.startsWith(line, "accession")) & !skipParse;
//					assert(valid); //Not true if concatenated
					
//					if(Tools.startsWith(line, "NZ_LM994619")){
//						boolean b=parseLine2(line, (byte)'\t');
//						assert(false) : b+", "+new String(line);
//					}
					
					if(valid){
						boolean b=parseLine2(line, (byte)'\t');
						if(b){linesValidT++;}
					}
				}
				line=bf.nextLine();
			}
			
			boolean closedError=bf.close();
			
//			if(USE_MAPS){
				for(int i=0; i<mapsT.length; i++){
					if(mapsT[i].size()>0){
						synchronized(maps[i]){
							maps[i].putAll(mapsT[i]);
						}
					}
					mapsT[i]=null;
				}
//			}
			if(USE_TABLES){
				long temp=table.flush();
			}
			
			success=!closedError;
		}
		
//		public boolean parseLineNumeric(final byte[] line, final byte delimiter){
//			int a=0, b=0;
//
//			long accession=0;
//			final int ncbi, gi;
//
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 0: "+new String(line);
//			for(int i=a; i<b; i++){
//				long c=line[i];
//				if(c=='.'){break;}
//				if(c<='9'){c=c-'0';}
//				else{c=c-'A'+10;}
//				accession=(accession*36)+c;
//			}
//			b++;
//			a=b;
//
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 1: "+new String(line);
//			//accession2=new String(line, a, b-a);
//			b++;
//			a=b;
//
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 2: "+new String(line);
//			ncbi=Parse.parseInt(line, a, b);
//			b++;
//			a=b;
//
////			while(b<line.length && line[b]!=delimiter){b++;}
////			assert(b>a) : "Missing field 3: "+new String(line);
//////			gi=Parse.parseInt(line, a, b);
////			b++;
////			a=b;
//
//			if(ncbi<1){return false;}
//
//			if(tree!=null){
//				if(ncbi>=tree.nodes.length){return false;}
//				TaxNode tn=tree.getNode(ncbi);
//				if(tn==null || tn.level==TaxTree.NO_RANK || tn.level==TaxTree.LIFE || tn.level==TaxTree.DOMAIN){return false;}
//				if(tn.pid>=tree.nodes.length){return false;}
//				tn=tree.getNode(tn.pid);
//				if(tn==null || tn.level==TaxTree.NO_RANK || tn.level==TaxTree.LIFE){return false;}
//			}
//			assert(accession>=0) : new String(line);
//			table.set(accession, ncbi);
//			return true;
//		}
		
		public boolean parseLine(final byte[] line, final byte delimiter){
			int a=0, b=0;
			
			String accession;
			final int ncbi, gi;
			
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing field 0: "+new String(line);
			accession=new String(line, a, b-a);
			final int dot=accession.indexOf('.');
			if(dot>=0){//Should never happen
//				System.err.println(accession);
//				assert(dot==accession.length()-2) : accession;
				accession=accession.substring(0, dot);
			}
//			if(STRIP_UNDERSCORE){
//				accession=accession.replaceAll("[_-]", "");
//			}
			if(lengthCountsT!=null){lengthCountsT[b-a]++;}
			if(symbolCountsT!=null){
				for(int i=a; i<b; i++){symbolCountsT[line[i]]++;}
			}
			final int underscore=accession.indexOf('_');
			if(underscore>=0){
				if(counts_underscoreT!=null){counts_underscoreT[b-a]++;}
				if(counts_underscore2T!=null && underscore==2){counts_underscore2T[b-a]++;}
			}
			b++;
			a=b;
			
			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 1: "+new String(line);
			assert(b>=a) : "Missing field 1: "+new String(line)+"\n"+a+", "+b;
			//accession2=new String(line, a, b-a);
			b++;
			a=b;
			
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing field 2: "+new String(line);
			ncbi=Parse.parseInt(line, a, b);
			b++;
			a=b;
			
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 3: "+new String(line);
////			gi=Parse.parseInt(line, a, b);
//			b++;
//			a=b;
			
			if(ncbi<1){return false;}
			
			if(tree!=null){
				if(ncbi>=tree.nodes.length){return false;}
				TaxNode tn=tree.getNode(ncbi);
				if(tn==null || tn.levelExtended==TaxTree.NO_RANK_E || tn.levelExtended==TaxTree.LIFE_E || tn.levelExtended==TaxTree.DOMAIN_E){return false;}
				if(tn.pid>=tree.nodes.length){return false;}
				tn=tree.getNode(tn.pid);
				if(tn==null || tn.levelExtended==TaxTree.NO_RANK_E || tn.levelExtended==TaxTree.LIFE_E){return false;}
			}
			
			if(accession.length()<13 && USE_TABLES){
				long number=hash(accession);
				assert(number>=0) : new String(line);
				table.set(number, ncbi);
				return true;
			}
			
			int way=accession.charAt(0);
			mapsT[way].put(accession, ncbi);
//			Integer old=mapsT[way].put(accession, ncbi);
//			assert(old==null || old==ncbi) : "'"+accession+"': "+old+" -> "+ncbi;
//			System.err.println("'"+accession+"': "+old+" -> "+ncbi);
//			assert(dot==-1) : "'"+accession+"': "+old+" -> "+ncbi;
			return true;
		}
		
		public boolean parseLine2(final byte[] line, final byte delimiter){
			int a=0, b=0;
			
			final int ncbi, gi;

			while(b<line.length && line[b]!=delimiter && line[b]!='.'){b++;}//parse unique part of accession
			final int dot=b;
			assert(b>a) : "Missing field 0: "+new String(line);
			while(b<line.length && line[b]!=delimiter){b++;}//skip the rest of the accession

			//System.err.println("Line: "+new String(line)+"\n"+Arrays.toString(line));
			//System.err.println("A: dot="+dot+", a="+a+", b="+b);
			
			{//Optional block
				if(lengthCountsT!=null){lengthCountsT[dot]++;}
				if(symbolCountsT!=null){
					for(int i=0; i<dot; i++){symbolCountsT[line[i]]++;}
				}
				if(counts_underscoreT!=null || counts_underscore2T!=null){
					int underscore=-1;
					for(int i=0; i<dot; i++){
						if(line[i]=='_'){
							underscore=i;
							break;
						}
					}
					if(underscore>=0){
						if(counts_underscoreT!=null){counts_underscoreT[dot]++;}
						if(counts_underscore2T!=null && underscore==2){counts_underscore2T[dot]++;}
					}
				}
			}
			b++;
			a=b;
			
			//System.err.println("B: a="+a+", b="+b);
			
			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 1: "+new String(line);
			assert(b>=a) : "Missing field 1: "+new String(line)+"\n"+a+", "+b;
			//accession2=new String(line, a, b-a);
			b++;
			a=b;
			
			//System.err.println("C: a="+a+", b="+b);
			
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing field 2: "+new String(line);
			ncbi=Parse.parseInt(line, a, b);
			//System.err.println("D: a="+a+", b="+b+", ncbi="+ncbi+", '"+(new String(line, a, b-a))+"'");
			b++;
			a=b;
			
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 3: "+new String(line);
////			gi=Parse.parseInt(line, a, b);
//			b++;
//			a=b;
			
			if(ncbi<1){return false;}
			//System.err.println("E: a="+a+", b="+b);
			if(skipHash){return false;}//123
			//System.err.println("F: a="+a+", b="+b);
			
			if(tree!=null){
				if(ncbi>=tree.nodes.length){return false;}
				//System.err.println("G");
				TaxNode tn=tree.getNode(ncbi);
				if(tn==null || /*tn.levelExtended==TaxTree.NO_RANK_E ||*/ tn.levelExtended==TaxTree.LIFE_E || tn.levelExtended==TaxTree.DOMAIN_E){return false;}
				//System.err.println("H: "+tn);
				if(tn.pid>=tree.nodes.length){return false;}
				//System.err.println("I: "+tn);
//				TaxNode parent=tree.getNode(tn.pid);
//				System.err.println("J: "+tn);
//				if(tn==null || tn.levelExtended==TaxTree.NO_RANK_E || tn.levelExtended==TaxTree.LIFE_E){return false;}
//				System.err.println("K");
			}
			
			if(USE_TABLES){
				if(AnalyzeAccession.codeMap!=null){
//					if(dot>AnalyzeAccession.longestPattern){return false;}
					final long number=AnalyzeAccession.digitize(line);
					if(number>=0){
						table.set(number, ncbi);
						return true;
					}
					assert(number==-1) : number+", "+new String(line);
				}else{
					if(dot<13){
						//				long number=hash(accession);
						final long number=hash(line, dot);
						assert(number>=0) : new String(line);
						table.set(number, ncbi);
						return true;
					}
				}
			}
			
			String accession=new String(line, 0, dot);
			int way=accession.charAt(0);
			mapsT[way].put(accession, ncbi);
//			Integer old=mapsT[way].put(accession, ncbi);
//			assert(old==null || old==ncbi) : "'"+accession+"': "+old+" -> "+ncbi;
//			System.err.println("'"+accession+"': "+old+" -> "+ncbi);
//			assert(dot==-1) : "'"+accession+"': "+old+" -> "+ncbi;
			return true;
		}
		
		private long linesProcessedT=0;
		private long linesValidT=0;
		private long bytesProcessedT=0;
		
		final FileFormat ff;
		HashMap<String, Integer>[] mapsT;
		HashBuffer table;
		boolean success=false;
		
		private long[] lengthCountsT=null;//new AtomicLongArray(20);
		private long[] symbolCountsT=null;//new AtomicLongArray(255);
		private long[] counts_underscoreT=null;//new AtomicLongArray(20);
		private long[] counts_underscore2T=null;//new AtomicLongArray(20);
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
//	private String out=null;
	
	static int maxPigzProcesses=12;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesValid=0;
	private long bytesProcessed=0;
	
	private AtomicLongArray lengthCounts=null;//new AtomicLongArray(20);
	private AtomicLongArray symbolCounts=null;//new AtomicLongArray(255);
	private AtomicLongArray counts_underscore=null;//new AtomicLongArray(20);
	private AtomicLongArray counts_underscore2=null;//new AtomicLongArray(20);
	
	/*--------------------------------------------------------------*/

	private final FileFormat ffin[];
//	private final FileFormat ffout;
	
	
	/*--------------------------------------------------------------*/
	
	public static boolean LOADED(){return LOADED;}
	
	private static boolean LOADED=false;
	private static HashMap<String, Integer>[] maps=null;
	private static KmerTableSet tables;
	public static TaxTree tree=null;
//	public static final boolean USE_MAPS=true;
	public static final boolean USE_TABLES=true;
//	public static boolean STRIP_UNDERSCORE=false;
	public static boolean skipParse=false;
	public static boolean skipHash=false;
	public static float prealloc=0;
	private static final long offset=-'A'+11;
	private static final long offsetLower=-'a'+11;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	
}
