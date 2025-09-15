package gff;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import aligner.Alignment;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.PGMTools;
import prok.ProkObject;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import tax.GiToTaxid;
import tax.TaxClient;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

public class CutGff implements Accumulator<CutGff.ProcessThread>  {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class.
		CutGff x=new CutGff(args);
		
		//Run the object.
		x.process(t);
		
		//Close the print stream if it was redirected.
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CutGff(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Shared.TRIM_READ_COMMENTS=Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;
		Read.VALIDATE_IN_CONSTRUCTOR=true;
		GffLine.parseAttributes=true;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;

			out=parser.out1;
		}
		
		if(alignRibo){
			//Load sequences
			ProkObject.loadConsensusSequenceFromFile(false, false);
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout=FileFormat.testOutput(out, FileFormat.PGM, null, true, overwrite, append, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		parser.overwrite=overwrite;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

//			outstream.println(arg+", "+a+", "+b);
			if(PGMTools.parseStatic(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna") || a.equals("ref")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(a.equals("gff") || a.equals("ingff") || a.equals("gffin")){
				assert(b!=null);
				Tools.addFiles(b, gffList);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("alignribo") || a.equals("align")){
				alignRibo=Parse.parseBoolean(b);
			}else if(a.equals("adjustendpoints")){
				adjustEndpoints=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("slop16s") || a.equalsIgnoreCase("16sslop") || a.equalsIgnoreCase("ssuslop")){
				ssuSlop=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("slop23s") || a.equalsIgnoreCase("23sslop") || a.equalsIgnoreCase("lsuslop")){
				lsuSlop=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxns") || a.equalsIgnoreCase("maxundefined")){
				maxNs=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxnrate") || a.equalsIgnoreCase("maxnfraction")){
				maxNFraction=Integer.parseInt(b);
			}else if(a.equals("invert")){
				invert=Parse.parseBoolean(b);
			}else if(a.equals("type") || a.equals("types")){
				types=b;
			}else if(a.equals("attributes") || a.equals("requiredattributes")){
				requiredAttributes=b.split(",");
			}else if(a.equals("banattributes") || a.equals("bannedattributes")){
				bannedAttributes=b.split(",");
			}else if(a.equals("banpartial")){
				banPartial=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("renameByTaxID")){
				renameByTaxID=Parse.parseBoolean(b);
			}else if(a.equals("taxmode")){
				if("accession".equalsIgnoreCase(b)){
					taxMode=ACCESSION_MODE;
				}else if("header".equalsIgnoreCase(b)){
					taxMode=HEADER_MODE;
				}else if("gi".equalsIgnoreCase(b)){
					taxMode=GI_MODE;
				}else if("taxid".equalsIgnoreCase(b)){
					taxMode=TAXID_MODE;
				}else{
					assert(false) : "Bad tax mode: "+b;
				}
			}else if(a.equals("requirepresent")){
				requirePresent=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("onePerFile")){
				onePerFile=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("pickBest") || a.equalsIgnoreCase("findBest") || a.equalsIgnoreCase("keepBest")){
				pickBest=Parse.parseBoolean(b);
			}
			
			else if(a.equals("minlen")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("maxlen")){
				maxLen=Integer.parseInt(b);
			}
			
			else if(ProkObject.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(arg.indexOf('=')<0 && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		ArrayList<String> banned=new ArrayList<String>();
		if(banPartial){banned.add("partial=true");}
		if(bannedAttributes!=null){
			for(String s : bannedAttributes){banned.add(s);}
		}
		bannedAttributes=banned.isEmpty() ? null : banned.toArray(new String[0]);
		
		if(gffList.isEmpty()){
			for(String s : fnaList){
				String prefix=ReadWrite.stripExtension(s);
				String gff=prefix+".gff";
				File f=new File(gff);
				if(!f.exists()){
					String gz=gff+".gz";
					f=new File(gz);
					assert(f.exists() && f.canRead()) : "Can't read file "+gff;
					gff=gz;
				}
				gffList.add(gff);
			}
		}
		assert(gffList.size()==fnaList.size()) : "Number of fna and gff files do not match: "+fnaList.size()+", "+gffList.size();
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		gffList=Tools.fixExtension(gffList);
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(gffList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		foo.add(out);
		if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Actual Code          ----------------*/
	/*--------------------------------------------------------------*/
	
	
	
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		if(Shared.threads()>2 && fnaList.size()>1){
			processMT(t);
		}else{
			processST(t);
		}
	}
	
	public void processST(Timer t){
		ByteStreamWriter bsw=(ffout==null ? null : new ByteStreamWriter(ffout));
		if(bsw!=null){bsw.start();}
		
		for(int i=0; i<fnaList.size(); i++){
			processFileST(fnaList.get(i), gffList.get(i), types, bsw);
		}
		
		if(bsw!=null){bsw.poisonAndWait();}
		t.stop();
		if(ffout!=null){outstream.println("Wrote "+out);}
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		if(alignRibo){
			outstream.println(Tools.number("Flipped:           ", flipped.get(), 8));
			outstream.println(Tools.number("Failed Alignment:  ", failed.get(), 8));
		}
	}
	
	public void processMT(Timer t){
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(ros);
		
		//Report timing and results
		t.stop();
		if(ffout!=null){outstream.println("Wrote "+out);}
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		if(alignRibo){
			outstream.println(Tools.number("Flipped:           ", flipped.get(), 8));
			outstream.println(Tools.number("Failed Alignment:  ", failed.get(), 8));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private boolean hasAttributes(GffLine gline){
		int len=gline.length();
		if(len<minLen || len>maxLen){return false;}
		if(hasAttributes(gline, bannedAttributes)){return false;}
		return requiredAttributes==null || hasAttributes(gline, requiredAttributes);
	}
	
	private static boolean hasAttributes(GffLine gline, String[] attributes){
		if(attributes==null){return false;}
		for(String s : attributes){
			if(gline.attributes.contains(s)){
				return true;
			}
		}
		return false;
	}
	
	private void processFileST(String fna, String gff, String types, ByteStreamWriter bsw){
		ArrayList<Read> reads=processFile(fna, gff, types);
		if(reads!=null){
			for(Read r : reads){
				if(bsw!=null){bsw.println(r);}
			}
		}
	}
	
	private ArrayList<Read> processFile(String fna, String gff, String types){
		ArrayList<GffLine> lines=GffLine.loadGffFile(gff, types, false);
		
		ArrayList<Read> list=ReadInputStream.toReads(fna, FileFormat.FA, -1);
		HashMap<String, Read> map=new HashMap<String, Read>();
		
		for(Read r : list){
			readsProcessed++;
			basesProcessed+=r.length();
			map.put(r.id, r);
		}
		
		if(renameByTaxID){//Note this must be AFTER adding to the hashmap.
			renameByTaxID(list);
		}
		
		ArrayList<Read> outList=processLines(lines, map, invert);
		
		if(invert){
			for(Read r : list){
				readsOut++;
				basesOut+=r.length();
			}
			return list;
		}else{
			if(outList!=null){
				for(Read r : outList){
					readsOut++;
					basesOut+=r.length();
				}
			}
			return outList;
		}
	}
	
	private void renameByTaxID(ArrayList<Read> list){
		ByteBuilder bb=new ByteBuilder();
		for(Read r : list){
			if(bb.length>0){bb.comma();}
			bb.append(taxMode==HEADER_MODE ? r.id : taxMode==ACCESSION_MODE ? parseAccession(r.id) : parseGi(r.id));
		}
		final int[] ids;
		if(taxMode==ACCESSION_MODE){
			ids=TaxClient.accessionToTaxidArray(bb.toString());
		}else if(taxMode==GI_MODE){
			ids=TaxClient.giToTaxidArray(bb.toString());
		}else if(taxMode==HEADER_MODE){
			ids=TaxClient.headerToTaxidArray(bb.toString());
		}else if(taxMode==TAXID_MODE){
			ids=parseTaxIds(list);
		}else{
			throw new RuntimeException("Bad mode: "+TAXID_MODE);
		}
		assert(ids.length==list.size()) : ids.length+", "+list.size();
		for(int i=0; i<ids.length; i++){
			Read r=list.get(i);
			int id=ids[i];
			
			if(r.id!=null && r.id.startsWith("tid|")){
				id=TaxTree.parseHeaderStatic(r.id);
				r.obj=id;
			}else {
				assert(id>=0 || !requirePresent) : "Can't find taxID for header: "+id+", "+r.name();

				r.obj=id;
				r.id="tid|"+id+"|"+id;
			}
		}
	}
	
	private int[] parseTaxIds(ArrayList<Read> list){
		int[] ids=new int[list.size()];
		for(int i=0; i<list.size(); i++){
			Read r=list.get(i);
			int x=GiToTaxid.parseTaxidNumber(r.id, '|');
			ids[i]=x;
		}
		return ids;
	}
	
	private String parseAccession(String id){
		int dot=id.indexOf('.');
		int space=id.indexOf(' ');
//		assert(dot>0 && space>0) : "Unexpected header format: "+id+"\nTry header mode instead of accession mode.";
		int limit=Tools.min((dot<0 ? id.length() : dot), (space<0 ? id.length() : space));
		return id.substring(0, limit);
	}
	
	private String parseGi(String id){
		assert(id.startsWith("gi|"));
		String[] split=id.split("|");
		return split[1];
	}
	
	private ArrayList<Read> processLines(ArrayList<GffLine> lines, HashMap<String, Read> map, boolean invertSelection){
		ArrayList<Read> list=null; 
		for(GffLine gline : lines){
			if(hasAttributes(gline)){
				Read scaf=map.get(gline.seqid);
				assert(scaf!=null) : "Can't find "+gline.seqid+" in "+map.keySet();
				
				boolean pass=true;
				Float identity=null;
				if(alignRibo && gline.inbounds(scaf.length())){
					int type=gline.prokType();
					identity=align(gline, scaf.bases, type);
					if(identity==null){pass=false;}
				}
				
				if(pass){
					final int start=gline.start-1;
					final int stop=gline.stop-1;

					if(invertSelection){
						byte[] bases=scaf.bases;
						for(int i=start; i<=stop; i++){
							if(i>=0 && i<bases.length){
								bases[i]='N';
							}
						}
					}else{
						if(start>=0 && stop<scaf.length()){
							String id=gline.attributes;
							if(renameByTaxID){
								id="tid|"+scaf.obj+"|"+id;
							}
							Read r=new Read(Arrays.copyOfRange(scaf.bases, start, stop+1), null, id, 1);
							r.obj=identity;
							
							assert(!r.containsLowercase()) : r.toFasta()+"\n"
							+ "validated="+r.validated()+", scaf.validated="+scaf.validated()+", tuc="+Read.TO_UPPER_CASE+", vic="+Read.VALIDATE_IN_CONSTRUCTOR;
							if(maxNs>=0 || maxNFraction>=0){
								long allowed=Tools.min(maxNs>=0 ? maxNs : r.length(), (long)(r.length()*(maxNFraction>=0 ? maxNFraction : 1)));
								if(r.countUndefined()>allowed){r=null;}
							}
							
							if(r!=null){
								if(gline.strand==1){r.reverseComplement();}
								if(list==null){list=new ArrayList<Read>(8);}
								list.add(r);
								if(onePerFile && !pickBest){break;}
							}
						}
					}
				}
			}
		}
		if(pickBest && list!=null && list.size()>1){
			Read best=null;
			float bestID=0;
			for(Read r : list){
				float identity=(r.obj==null ? 0.001f : (Float)r.obj);
				if(identity>bestID){
					bestID=identity;
					best=r;
				}
			}
			assert(best!=null);
			list.clear();
			list.add(best);
		}
		return list;
	}
	
	private Float align(GffLine gline, byte[] scaf, int type){
		Read[] consensusReads=ProkObject.consensusReads(type);
		if(consensusReads==null || consensusReads.length==0){
			assert(false) : type+"\n"+gline.toString();
			return null;
		}
		byte[] universal=consensusReads[0].bases;
		float minIdentity=ProkObject.minID(type)*ID_MULT;
		if(universal==null){assert(false); return 1F;}
		
		int start=gline.start-1;
		int stop=gline.stop-1;
		assert(start<=stop) : start+", "+stop+", "+scaf.length;
		assert(start>=0 && start<scaf.length) : start+", "+stop+", "+scaf.length;
		assert(stop>=0 && stop<scaf.length) : start+", "+stop+", "+scaf.length;
//		final int a=Tools.max(0, start-(adjustEndpoints ? 100 : 20));
//		final int b=Tools.min(scaf.length-1, stop+(adjustEndpoints ? 100 : 20));
		final int a=Tools.max(0, start);
		final int b=Tools.min(scaf.length-1, stop);
		
		byte[] ref=Arrays.copyOfRange(scaf, a, b+1);
		Read r=new Read(ref, null, 0);
		if(gline.strand==GffLine.MINUS){r.reverseComplement();}

		Alignment plus=new Alignment(r);
		plus.align(universal);
//		if(plus.id>=minIdentity){return true;}

		r.reverseComplement();
		Alignment minus=new Alignment(r);
		minus.align(universal);
		
		Alignment best=null;
		if(plus.id>=minus.id){
			best=plus;
//			System.err.println("Kept:    "+plus.id+" \t"+minus.id);
		}else{
			best=minus;
			if(minus.id>=minIdentity){
				if(verbose) {System.err.println("Flipped: "+plus.id+" \t"+minus.id+"");}
				flipped.incrementAndGet();
				gline.strand=Shared.MINUS;
			}
		}
		if(best.id>=minIdentity){
			return best.id;
		}else{
			if(verbose) {System.err.println("Failed alignment: "+plus.id+" \t"+minus.id);}
			failed.incrementAndGet();
			return null;
		}
//		
//		int[] coords=KillSwitch.allocInt1D(2);
//		float id1=align(universal, scaf, a, b, minIdentity, coords);
//		final int rstart=coords[0], rstop=coords[1];
//		//				assert(false) : rstart+", "+rstop+", "+(rstop-rstart+1)+", "+start+", "+stop;
//		if(id1<minIdentity){
//			System.err.println("Low identity: "+String.format("%.2s", 100*id1));
//			return false;
//		}
//		if(adjustEndpoints){
//			int slop=(flag==4 ? AnalyzeGenes.lsuSlop : AnalyzeGenes.ssuSlop);
//			if(Tools.absdif(start, rstart)>slop){
//				//						System.err.println("rstart:\t"+start+" -> "+rstart);
//				start=rstart;
//			}
//			if(Tools.absdif(stop, rstop)>slop){
//				//						System.err.println("rstop: \t"+stop+" -> "+rstop);
//				stop=rstop;
//			}
//		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	private ConcurrentReadOutputStream makeCros(){
		if(ffout==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/** Spawn process threads */
	private void spawnThreads(ConcurrentReadOutputStream ros){
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Tools.min(Shared.threads(), fnaList.size());
		
		//Controls access to input files
		AtomicInteger atom=new AtomicInteger(0);
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(atom, ros, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		readsOut+=pt.readsOutT;
		basesOut+=pt.basesOutT;
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final AtomicInteger atom_, ConcurrentReadOutputStream ros_, final int tid_){
			atom=atom_;
			ros=ros_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			for(int fnum=atom.getAndIncrement(), lim=fnaList.size(); fnum<lim; fnum=atom.getAndIncrement()) {
				String fna=fnaList.get(fnum);
				String gff=gffList.get(fnum);
				//Process the reads
				ArrayList<Read> list=processFileT(fna, gff, types);
				if(ros!=null){
					if(list==null){list=dummy;}
					ros.add(list, fnum);
				}
			}
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		//Duplicated
		private ArrayList<Read> processFileT(String fna, String gff, String types){
			ArrayList<GffLine> lines=GffLine.loadGffFile(gff, types, false);
			
			ArrayList<Read> list=ReadInputStream.toReads(fna, FileFormat.FA, -1);
			HashMap<String, Read> map=new HashMap<String, Read>();
			
			for(Read r : list){
				readsProcessedT++;
				basesProcessedT+=r.length();
				map.put(r.id, r);
			}
			
			if(renameByTaxID){//Note this must be AFTER adding to the hashmap.
				renameByTaxID(list);
			}
//			assert(false) : renameByTaxID+", "+list.size();
			
			ArrayList<Read> outList=processLines(lines, map, invert);
			
			if(invert){
				for(Read r : list){
					readsOutT++;
					basesOutT+=r.length();
				}
				return list;
			}else{
				if(outList!=null){
					for(Read r : outList){
						readsOutT++;
						basesOutT+=r.length();
					}
				}
				return outList;
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
		/** Next file ID */
		final AtomicInteger atom;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> gffList=new ArrayList<String>();
	private String out=null;
	private String types="CDS";
	private boolean invert=false;
	private boolean banPartial=true;
	private int minLen=1;
	private int maxLen=Integer.MAX_VALUE;

	private String[] requiredAttributes;
	private String[] bannedAttributes;
	
	/*--------------------------------------------------------------*/
	
	private long bytesOut=0;
	private boolean renameByTaxID=false;
	private int taxMode=ACCESSION_MODE;
	private boolean requirePresent=false;
	private boolean alignRibo=false;
	private boolean adjustEndpoints=false;
	private boolean onePerFile=false;
	private boolean pickBest=false;
	private int ssuSlop=999;
	private int lsuSlop=999;
	
	private float ID_MULT=0.96f;
	
	private int maxNs=-1;
	private double maxNFraction=-1;
	
	private static int ACCESSION_MODE=0, GI_MODE=1, HEADER_MODE=2, TAXID_MODE=3;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	protected AtomicLong flipped=new AtomicLong(0);
	protected AtomicLong failed=new AtomicLong(0);

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffout;
	private final ArrayList<Read> dummy=new ArrayList<Read>();
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	public boolean ordered=true;
	private boolean overwrite=true;
	private boolean append=false;
	
}
