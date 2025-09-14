package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.LongListSet;
import structures.LongListSet.LongListSetIterator;
import tracker.ReadStats;

/**
 * Generates mutants of kmers.
 * 
 * @author Brian Bushnell
 * @date January 8, 2021
 *
 */
public class KExpand {
	
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
		
		//Create an instance of this class
		KExpand x=new KExpand(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KExpand(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		Shared.capBuffers(4); //Only for singlethreaded programs
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			extin=parser.extin;

			out1=parser.out1;
			extout=parser.extout;
			amino=Shared.AMINO_IN;
		}

		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program
		
//		k2=k-1;
//		shift=bitsPerBase*k;
//		shift2=shift-bitsPerBase;
//		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		
		{//set some constants
			k2=k-1;
			bitsPerBase=(amino ? 5 : 2);
			maxSymbol=(amino ? 20 : 3);
			symbols=maxSymbol+1;
			symbolArrayLen=(64+bitsPerBase-1)/bitsPerBase;
			symbolSpace=(1<<bitsPerBase);
			symbolMask=symbolSpace-1;
			
			symbolToNumber=AminoAcid.symbolToNumber(amino);
			symbolToNumber0=AminoAcid.symbolToNumber0(amino);
			symbolToComplementNumber0=AminoAcid.symbolToComplementNumber0(amino);
			
			clearMasks=new long[symbolArrayLen];
			leftMasks=new long[symbolArrayLen];
			rightMasks=new long[symbolArrayLen];
//			lengthMasks=new long[symbolArrayLen];
			setMasks=new long[symbols][symbolArrayLen];
			for(int i=0; i<symbolArrayLen; i++){
				clearMasks[i]=~(symbolMask<<(bitsPerBase*i));
				leftMasks[i]=((-1L)<<(bitsPerBase*i));
				rightMasks[i]=~((-1L)<<(bitsPerBase*i));
//				lengthMasks[i]=((1L)<<(bitsPerBase*i));
				for(long j=0; j<symbols; j++){
					setMasks[(int)j][i]=(j<<(bitsPerBase*i));
				}
			}
			
			minlen=k-1;
			minlen2=k;
			shift=bitsPerBase*k;
			shift2=shift-bitsPerBase;
			mask=(shift>63 ? -1L : ~((-1L)<<shift));
//			kmask=lengthMasks[k];
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		parser.overwrite=true;
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Parse.parseIntKMG(b);
				assert(k>0 && k<=31);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("subdist") || a.equals("sdist") || a.equals("hdist") || a.equals("hammingdistance")){
				subDist=Parse.parseIntKMG(b);
				assert(subDist>=0);
			}else if(a.equals("deldist") || a.equals("ddist") || a.equals("deletiondistance")){
				delDist=Parse.parseIntKMG(b);
				assert(delDist>=0 && delDist<=3);
			}else if(a.equals("insdist") || a.equals("idist") || a.equals("insertiondistance")){
				insDist=Parse.parseIntKMG(b);
				assert(insDist>=0);
			}else if(a.equals("editdist") || a.equals("edist") || a.equals("edits") || a.equals("editdistance")){
				int x=Parse.parseIntKMG(b);
				assert(x>=0 && x<=3);
				editDist=x;
				subDist=Tools.max(subDist, x);
			}else if(a.equals("maxedits") || a.equals("emax")){
				maxEdits=Parse.parseIntKMG(b);
				assert(maxEdits>=0);
			}else if(a.equals("maxsubs") || a.equals("smax")){
				maxSubs=Parse.parseIntKMG(b);
				assert(maxSubs>=0);
			}else if(a.equals("maxdels") || a.equals("dmax")){
				maxDels=Parse.parseIntKMG(b);
				assert(maxDels>=0);
			}else if(a.equals("maxinss") || a.equals("imax")){
				maxInss=Parse.parseIntKMG(b);
				assert(maxInss>=0);
			}else if(a.equals("speed")){
				int x=Parse.parseIntKMG(b);
				assert(x>=0 && x<=16) : "Speed must be 0-16";
				speed=x;
			}else if(a.equals("skip")){
				int x=Parse.parseIntKMG(b);
				assert(x>=0) : "Skip must be >=0";
				skip=x;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Make sure interleaving agrees with number of input and output files */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the read stream
		processInner(cris);
		errorState|=ReadWrite.closeStreams(cris);

		if(verbose){outstream.println("Finished reading input.");}
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		if(ros!=null){
			dumpKmers(ros);
			if(verbose){outstream.println("Finished processing output.");}
			
			errorState|=ReadWrite.closeStream(ros);
		}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.things("Kmers Processed", kmersProcessed, 8));
		outstream.println(Tools.things("Kmers Added", kmersAdded, 8));
		
		outstream.println();
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		outstream.println(Tools.things("Kmers Out", kmersOut, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Dumping           ----------------*/
	/*--------------------------------------------------------------*/
	
	private void dumpKmers(ConcurrentReadOutputStream ros) {
		
		set.sort();
		set.shrinkToUnique();
		LongListSetIterator it=set.iterator();
		
		long id=1;
		ArrayList<Read> list=new ArrayList<Read>(200);
		ByteBuilder bb=new ByteBuilder();
		while(it.hasMore()){
			bb.clear();
			final long kmer=it.next();
			bb.appendKmer(kmer, k);
			Read r=new Read(bb.toBytes(), null, id);
			id++;
			kmersOut++;
			readsOut++;
			basesOut+=r.length();
			list.add(r);
			if(list.size()>=200){
				ros.add(list, 0);
				list=new ArrayList<Read>(200);
			}
		}
		if(list.size()>0){ros.add(list, 0);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	private ConcurrentReadOutputStream makeCros(){
		if(ffout1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=4;

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired());
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln, cris);

				//Fetch a new list
				ln=cris.nextList();
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
		
	}
	
	/**
	 * Process a list of Reads.
	 * @param ln The list.
	 * @param cris Read Input Stream
	 * @param ros Read Output Stream for reads that will be retained
	 */
	void processList(ListNum<Read> ln, final ConcurrentReadInputStream cris){

		//Grab the actual read list from the ListNum
		final ArrayList<Read> reads=ln.list;
		
		long added=0;
		
		//Loop through each read in the list
		for(int idx=0; idx<reads.size(); idx++){
			final Read r1=reads.get(idx);
			final Read r2=r1.mate;
			
			//Validate reads in worker threads
			if(!r1.validated()){r1.validate(true);}
			if(r2!=null && !r2.validated()){r2.validate(true);}

			//Track the initial length for statistics
			final int initialLength1=r1.length();
			final int initialLength2=r1.mateLength();

			//Increment counters
			readsProcessed+=r1.pairCount();
			basesProcessed+=initialLength1+initialLength2;
			kmersAdded+=processReadPair(r1, r2);
		}

		//Notify the input stream that the list was used
		cris.returnList(ln);
//		if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
	}
	
	
	/**
	 * Process a single read pair.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return Number of kmers added.
	 */
	long processReadPair(final Read r1, final Read r2){
		long x=addToMap(r1, skip);
		if(r2!=null){
			x+=addToMap(r2, skip);
		}
		return x;
	}
	
	private long addToMap(Read r, int skip){
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		long added=0;
		int len=0;
		
//		if(bases!=null){
//			readsProcessed++;
//			basesProcessed+=bases.length;
//		}
		if(bases==null || bases.length<k){return 0;}
		
		if(skip>1){ //Process while skipping some kmers
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=symbolToNumber0[b];
				long x2=symbolToComplementNumber0[b];
				kmer=((kmer<<bitsPerBase)|x)&mask;
				rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
				if(isFullyDefined(b)){len++;}else{len=0; rkmer=0;}
				if(verbose){outstream.println("Scanning1 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					kmersProcessed++;
					if(len%skip==0){
						final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
						final long extraBase2=(i>=bases.length-2 ? -1 : symbolToNumber[bases[i+2]]);
						final long extraBase3=(i>=bases.length-3 ? -1 : symbolToNumber[bases[i+3]]);
						added+=addAndMutate(kmer, rkmer, k, extraBase, extraBase2, extraBase3);
					}
				}
			}
		}else{ //Process all kmers
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=symbolToNumber0[b];
				final long x2=symbolToComplementNumber0[b];
//				assert(x!=x2) : x+", "+x2+", "+Character.toString((char)b)+"\n"+Arrays.toString(symbolToNumber0)+"\n"+Arrays.toString(symbolToComplementNumber);
				kmer=((kmer<<bitsPerBase)|x)&mask;
				//10000, 1111111111, 16, 16, 2, 10, 8
				rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
				if(isFullyDefined(b)){len++;}else{len=0; rkmer=0;}
				if(verbose){
//					if(verbose){
//						String fwd=new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k));
//						String rev=AminoAcid.reverseComplementBases(fwd);
//						String fwd2=kmerToString(kmer, Tools.min(len, k));
//						outstream.println("fwd="+fwd+", fwd2="+fwd2+", rev="+rev+", kmer="+kmer+", rkmer="+rkmer);
//						outstream.println("b="+(char)b+", x="+x+", x2="+x2+", bitsPerBase="+bitsPerBase+", shift2="+shift2);
//						if(!amino){
//							assert(AminoAcid.stringToKmer(fwd)==kmer) : fwd+", "+AminoAcid.stringToKmer(fwd)+", "+kmer+", "+len;
//							if(len>=k){
//								assert(rcomp(kmer, Tools.min(len, k))==rkmer);
//								assert(rcomp(rkmer, Tools.min(len, k))==kmer);
//								assert(AminoAcid.kmerToString(kmer, Tools.min(len, k)).equals(fwd));
//								assert(AminoAcid.kmerToString(rkmer, Tools.min(len, k)).equals(rev)) : AminoAcid.kmerToString(rkmer, Tools.min(len, k))+" != "+rev+" (rkmer)";
//							}
//							assert(fwd.equalsIgnoreCase(fwd2)) : fwd+", "+fwd2; //may be unsafe
//						}
//						outstream.println("Scanning6 i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+", bases="+fwd+", rbases="+rev);
//					}
				}
				if(len>=k){
//					assert(kmer==rcomp(rkmer, k)) : Long.toBinaryString(kmer)+", "+Long.toBinaryString(rkmer)+", "+Long.toBinaryString(mask)+", x="+x+", x2="+x2+", bits="+bitsPerBase+", s="+shift+", s2="+shift2+", b="+Character.toString((char)b);
					kmersProcessed++;
					final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
					final long extraBase2=(i>=bases.length-2 ? -1 : symbolToNumber[bases[i+2]]);
					final long extraBase3=(i>=bases.length-3 ? -1 : symbolToNumber[bases[i+3]]);
					final long atm=addAndMutate(kmer, rkmer, k, extraBase, extraBase2, extraBase3);
					added+=atm;
				}
			}
		}
		return added;
	}
	

	
	/**
	 * Adds this kmer to the table, including any mutations implied by subDist, etc.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param len Kmer length
	 * @param subDist Number of substitutions to allow
	 * @param delDist Number of deletions to allow
	 * @param insDist Number of insertions to allow
	 * @param extraBase Base added to end in case of deletions
	 * @param extraBase2 Base added to end in case of deletions
	 * @param extraBase3 Base added to end in case of deletions
	 * @return Number of kmers stored
	 */
	private long addAndMutate(final long kmer, final long rkmer, final int len, final long extraBase, final long extraBase2, final long extraBase3){
		if(editDist>0){//Edit mode
			return mutateE(kmer, rkmer, len, editDist, maxSubs, maxDels, maxInss, extraBase, extraBase2, extraBase3);
		}else if(subDist>0 || delDist>0 || insDist>0){//Sub Del Ins mode
			return mutateE(kmer, rkmer, len, maxEdits, subDist, delDist, insDist, extraBase, extraBase2, extraBase3);
		}else{
			addKey(kmer, rkmer, len);
		}
		return 1;
	}
	
	
	/**
	 * Adds this kmer to the set.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param len Kmer length
	 * @return Number of kmers added
	 */
	private long addKey(final long kmer, final long rkmer, final int len){
		if(verbose){outstream.println("addToMap_A; len="+len);}

		final long key=toValue(kmer, rkmer);
		if(verbose){outstream.println("toValue ("+kmerToString(kmer, len)+", "+kmerToString(rkmer, len)+") = "+kmerToString(key, len)+" = "+key);}
		if(failsSpeed(key)){return 0;}
		
		if(verbose){outstream.println("addToMap_B: "+kmerToString(key, len)+" ("+key+")");}
		long added=0;
		set.add(key);
		added++;
		if(verbose){outstream.println("addToMap added "+added+" keys.");}
		return 1;
	}
	
	//*** Note!  This is functionally identical to mutateE, 
	//*** only the default params are different
	@Deprecated
	/**
	 * Mutate and store this kmer through 'dist' recursions.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param len Kmer length
	 * @param subDist Number of substitutions to allow
	 * @param delDist Number of deletions to allow
	 * @param insDist Number of insertions to allow
	 * @param extraBase Base added to end in case of deletions
	 * @param extraBase2 Base added to end in case of deletions
	 * @param extraBase3 Base added to end in case of deletions
	 * @return Number of kmers stored
	 */
	private long mutateSDI(final long kmer, final long rkmer, final int len, final int maxEdits,
			final int subDist, final int delDist, final int insDist, 
			final long extraBase, final long extraBase2, final long extraBase3){
		long added=1;
		final int maxEdits2=maxEdits-1;
		
		addKey(kmer, rkmer, len);
		if(maxEdits<1){return 1;}
		
		if(verbose){outstream.println("mutateSDI; len="+len+"; kmer="+kmer+"; rkmer="+rkmer);}
		
		if(subDist>0){
			//Sub
			for(int j=0; j<symbols; j++){
				for(int i=0; i<len; i++){
					final long temp=(kmer&clearMasks[i])|setMasks[j][i];
					if(temp!=kmer){
						long rtemp=rcomp(temp, len);
						added+=mutateSDI(temp, rtemp, len, maxEdits2, 
								subDist-1, delDist, insDist, extraBase, extraBase2, extraBase3);
					}
				}
			}
		}
		
		if(delDist>0 && extraBase>=0 && extraBase<=maxSymbol){
			//Del
			for(int i=1; i<len; i++){
				final long temp=(kmer&leftMasks[i])|((kmer<<bitsPerBase)&rightMasks[i])|extraBase;
				if(temp!=kmer){
					long rtemp=rcomp(temp, len);
					added+=mutateSDI(temp, rtemp, len, maxEdits2, 
							subDist, delDist-1, insDist, extraBase2, extraBase3, -1);
				}
			}
		}
		
		if(insDist>0){
			//Ins
			final long eb0=kmer&symbolMask;
			for(int i=0; i<len; i++){
				final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>bitsPerBase);
				for(int j=0; j<symbols; j++){
					final long temp=temp0|setMasks[j][i];
					if(temp!=kmer){
						long rtemp=rcomp(temp, len);
						added+=mutateSDI(temp, rtemp, len, maxEdits2, 
								subDist, delDist, insDist-1, eb0, extraBase, extraBase2);
					}
				}
			}
		}
		
		return added;
	}
	
	/**
	 * Mutate and store this kmer through 'dist' recursions.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param len Kmer length
	 * @param editDist Number of edits to allow
	 * @param extraBase Base added to end in case of deletions
	 * @param extraBase2 Base added to end in case of deletions
	 * @param extraBase3 Base added to end in case of deletions
	 * @return Number of kmers stored
	 */
	private long mutateE(final long kmer, final long rkmer, final int len, final int editDist, 
			final int maxSubs, final int maxDels, final int maxInss, 
			final long extraBase, final long extraBase2, final long extraBase3){
		long added=0;
		
		addKey(kmer, rkmer, len);
		if(editDist<1){return 1;}
		
		if(verbose){outstream.println("mutateE; len="+len+"; kmer="+kmer+"; rkmer="+rkmer);}
		
		final int edist2=editDist-1;
		if(maxSubs>0){
			//Sub
			for(int j=0; j<symbols; j++){
				for(int i=0; i<len; i++){
					final long temp=(kmer&clearMasks[i])|setMasks[j][i];
					if(temp!=kmer){
						long rtemp=rcomp(temp, len);
						added+=mutateE(temp, rtemp, len, edist2, 
								maxSubs-1, maxDels, maxInss, extraBase, extraBase2, extraBase3);
					}
				}
			}
		}

		if(maxDels>0 && extraBase>=0 && extraBase<=maxSymbol){
			//Del
			for(int i=1; i<len; i++){
				final long temp=(kmer&leftMasks[i])|((kmer<<bitsPerBase)&rightMasks[i])|extraBase;
				if(temp!=kmer){
					long rtemp=rcomp(temp, len);
					added+=mutateE(temp, rtemp, len, edist2, 
							maxSubs, maxDels-1, maxInss, extraBase2, extraBase3, -1);
				}
			}
		}
		
		if(maxInss>0){
			//Ins
			final long eb0=kmer&symbolMask;
			for(int i=0; i<len; i++){
				final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>bitsPerBase);
				for(int j=0; j<symbols; j++){
					final long temp=temp0|setMasks[j][i];
					if(temp!=kmer){
						long rtemp=rcomp(temp, len);
						added+=mutateE(temp, rtemp, len, edist2,
								maxSubs, maxDels, maxInss-1, eb0, extraBase, extraBase2);
					}
				}
			}
		}
		
		return added;
	}
	
	/**
	 * Transforms a kmer into a canonical value.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @return Canonical value
	 */
	final long toValue(long kmer, long rkmer){
		if(verbose){outstream.println("toValue("+AminoAcid.kmerToString(kmer, k)+", "+AminoAcid.kmerToString(rkmer, k)+")");}
		final long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		if(verbose){outstream.println("value="+AminoAcid.kmerToString(value, k)+" = "+value);}
		return value;
	}
	
	final long rcomp(long kmer, int len){
		return amino ? kmer : AminoAcid.reverseComplementBinaryFast(kmer, len);
	}
	
	final boolean passesSpeed(long key){
		return speed<1 || ((key&Long.MAX_VALUE)%17)>=speed;
	}
	
	final boolean failsSpeed(long key){
		return speed>0 && ((key&Long.MAX_VALUE)%17)<speed;
	}
	
	/** Returns true if the symbol is not degenerate (e.g., 'N') for the alphabet in use. */
	final boolean isFullyDefined(byte symbol){
		return symbol>=0 && symbolToNumber[symbol]>=0;
	}
	
	/** For verbose / debugging output */
	final String kmerToString(long kmer, int k){
		return amino ? AminoAcid.kmerToStringAA(kmer, k) : AminoAcid.kmerToString(kmer, k);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;
	
	private boolean rcomp=true;
	private int k=31;
	/** k-1; used in some expressions */
	final int k2;

	private int speed=0;
	private int skip=0;

	private int editDist=0;
	private int maxSubs=99;
	private int maxDels=99;
	private int maxInss=99;

	private int maxEdits=99;
	private int subDist=0;
	private int delDist=0;
	private int insDist=0;
	
	private LongListSet set=new LongListSet();
	
	/*--------------------------------------------------------------*/
	/*-----------        Symbol-Specific Constants        ----------*/
	/*--------------------------------------------------------------*/

	/** True for amino acid data, false for nucleotide data */
	final boolean amino;
//	final int maxSupportedK;
	final int bitsPerBase;
	final int maxSymbol;
	final int symbols;
	final int symbolArrayLen;
	final int symbolSpace;
	final long symbolMask;
	
	final int minlen;
	final int minlen2;
	final int shift;
	final int shift2;
	final long mask;
	
	/** x&clearMasks[i] will clear base i */
	final long[] clearMasks;
	/** x|setMasks[i][j] will set base i to j */
	final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	final long[] rightMasks;
//	/** x|kMasks[i] will set the bit to the left of the leftmost base */
//	final long[] lengthMasks;
	
	private final byte[] symbolToNumber0;
	private final byte[] symbolToComplementNumber0;
	private final byte[] symbolToNumber;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of kmers processed */
	protected long kmersProcessed=0;
	/** Number of kmers added (includes mutants) */
	protected long kmersAdded=0;

	/** Number of reads output */
	protected long readsOut=0;
	/** Number of bases output */
	protected long basesOut=0;
	/** Number of kmers output */
	protected long kmersOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}
