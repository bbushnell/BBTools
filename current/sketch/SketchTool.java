package sketch;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerTableSet;
import parse.Parse;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.LongList;
import structures.StringNum;

/**
 * @author Brian Bushnell
 * @date June 28, 2016
 *
 */
public final class SketchTool extends SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a SketchTool with parameters from a DisplayParams object.
	 * @param size_ Target sketch size (number of k-mers to retain)
	 * @param params Display parameters containing configuration values
	 */
	public SketchTool(int size_, DisplayParams params){
		this(size_, params.minKeyOccuranceCount, params.trackCounts(), params.mergePairs, SketchObject.rcomp);
	}
	
	/**
	 * Main constructor for SketchTool with explicit parameters.
	 *
	 * @param size_ Target sketch size (number of k-mers to retain)
	 * @param minKeyOccuranceCount_ Minimum k-mer count threshold for inclusion
	 * @param trackCounts_ Whether to track k-mer occurrence counts
	 * @param mergePairs_ Whether to merge read pairs before sketching
	 * @param rcomp_ Whether to include reverse complement k-mers
	 */
	public SketchTool(int size_, int minKeyOccuranceCount_, boolean trackCounts_, 
			boolean mergePairs_, boolean rcomp_){
		stTargetSketchSize=size_;
		minKeyOccuranceCount=minKeyOccuranceCount_;
		trackCounts=trackCounts_;
		mergePairs=mergePairs_;
		rcomp=rcomp_;
		
		assert(!aminoOrTranslate() || !rcomp) : "rcomp should be false in amino mode.";
		assert(!aminoOrTranslate() || (k*AminoAcid.AMINO_SHIFT<64)) : "Protein sketches require 1 <= K <= "+(63/AminoAcid.AMINO_SHIFT)+".";
		assert(k>0 && k<=32) : "Sketches require 1 <= K <= 32."; //123
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Converts k-mer table set to a sketch using single or multithreaded approach.
	 * @param tables K-mer table set containing counted k-mers
	 * @param multithreaded Whether to use multithreaded processing
	 * @return Generated sketch from the k-mer tables
	 */
	public Sketch toSketch(KmerTableSet tables, boolean multithreaded){
		final int threads=(multithreaded ? Tools.mid(1, Shared.threads(), tables.ways()) : 1);
		return (threads<2 ? toSketch_ST(tables) : toSketch_MT(tables, threads));
	}
	
	/**
	 * Single-threaded sketch generation from k-mer tables.
	 * Creates either a fixed-size heap or unlimited list based on target size.
	 * @param tables K-mer table set to process
	 * @return Generated sketch
	 */
	private Sketch toSketch_ST(KmerTableSet tables){
		SketchHeap heap=(stTargetSketchSize>0 ? new SketchHeap(stTargetSketchSize, minKeyOccuranceCount, trackCounts) : null);
		LongList list=new LongList();
		
		KmerTableSet kts=(KmerTableSet)tables;
		for(int tnum=0; tnum<kts.ways; tnum++){
			HashArray1D table=kts.getTable(tnum);
			if(stTargetSketchSize>0){
				toHeap(table, heap);
			}else{
				toList(table, list);
			}
		}
		return stTargetSketchSize>0 ? new Sketch(heap, false, trackCounts, null) : toSketch(list);//TODO:  Could add counts here
	}
	
	/**
	 * Multithreaded sketch generation from k-mer tables.
	 * Distributes tables across worker threads for parallel processing.
	 *
	 * @param tables K-mer table set to process
	 * @param threads Number of worker threads to use
	 * @return Combined sketch from all threads
	 */
	private Sketch toSketch_MT(KmerTableSet tables, final int threads){
		ArrayList<SketchThread> alst=new ArrayList<SketchThread>(threads);
		AtomicInteger ai=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			alst.add(new SketchThread(ai, tables));
		}

		//Start the threads
		for(SketchThread pt : alst){
			pt.start();
		}

		ArrayList<SketchHeap> heaps=new ArrayList<SketchHeap>(threads);
		LongList list=new LongList();
		
		for(SketchThread pt : alst){

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
			
			if(stTargetSketchSize>=0){
				if(pt.heap!=null && pt.heap.size()>0){
					heaps.add(pt.heap);
				}
			}else{
				if(pt.list!=null){list.append(pt.list);}
				pt.list=null;
			}
		}
		alst.clear();
		
		return stTargetSketchSize>=0 ? toSketch(heaps, true) : toSketch(list);//TODO:  Could add counts here
	}
	
	/**
	 * Converts a single hash table to a sketch heap.
	 * Processes both main array and overflow forest structures.
	 *
	 * @param table Hash table containing k-mers and counts
	 * @param heap Target heap to accumulate k-mers (or null to create)
	 * @return The heap with added k-mers
	 */
	public SketchHeap toHeap(HashArray1D table, SketchHeap heap){
//		if(heap==null){heap=new LongHeap(size, true);}
		long[] kmers=table.array();
		int[] counts=table.values();
		for(int i=0; i<table.arrayLength(); i++){
			int count=counts[i];
			if(count>=minKeyOccuranceCount){
				heap.genomeSizeKmers++;
				long kmer=kmers[i];
				long hashcode=hash(kmer);
				if(hashcode>=minHashValue){
					heap.add(hashcode);
				}
			}
		}
		HashForest forest=table.victims();
		if(forest!=null){
			for(KmerNode kn : forest.array()){
				if(kn!=null){addRecursive(heap, kn);}
			}
		}
		return heap;
	}
	
	/**
	 * Converts a hash table to an unlimited list of k-mer hashes.
	 * Used when no size limit is specified for the sketch.
	 *
	 * @param table Hash table containing k-mers and counts
	 * @param list Target list to accumulate k-mers
	 * @return The list with added k-mers
	 */
	public LongList toList(HashArray1D table, LongList list){
//		if(heap==null){heap=new LongHeap(size, true);}
		long[] kmers=table.array();
		int[] counts=table.values();
		for(int i=0; i<table.arrayLength(); i++){
			int count=counts[i];
			if(count>=minKeyOccuranceCount){
				long kmer=kmers[i];
				long hashcode=hash(kmer);
				if(hashcode>=minHashValue){
					list.add(hashcode);
				}
			}
		}
		HashForest forest=table.victims();
		if(forest!=null){
			for(KmerNode kn : forest.array()){
				if(kn!=null){addRecursive(list, kn);}
			}
		}
		return list;
	}
	
//	public long[] toSketchArray(ArrayList<LongHeap> heaps){
//		if(heaps.size()==1){return toSketchArray(heaps.get(0));}
//		LongList list=new LongList(size);
//		for(LongHeap heap : heaps){
//			while(heap.size()>0){list.add(Long.MAX_VALUE-heap.poll());}
//		}
//		list.sort();
//		list.shrinkToUnique();
//		list.size=Tools.min(size, list.size);
//		return list.toArray();
//	}
	
	/**
	 * Combines multiple sketch heaps into a single sketch.
	 * Merges heaps by adding them together sequentially.
	 *
	 * @param heaps List of sketch heaps to combine
	 * @param allowZeroSizeSketch Whether to allow creation of empty sketches
	 * @return Combined sketch or null if no valid heaps
	 */
	public Sketch toSketch(ArrayList<SketchHeap> heaps, boolean allowZeroSizeSketch){
		if(heaps==null || heaps.isEmpty()){
			if(allowZeroSizeSketch){
				return new Sketch(new long[0], null, null, null, null, null);
			}else{
				return null;
			}
		}
		SketchHeap a=heaps.get(0);
		for(int i=1; i<heaps.size(); i++){
			SketchHeap b=heaps.get(i);
			a.add(b);
		}
		if(verbose2){System.err.println("Creating a sketch of size "+a.size()+".");}
		return new Sketch(a, false, trackCounts, null);
	}
	
	/**
	 * Converts a list of k-mer hashes to a sketch.
	 * Sorts, deduplicates, and transforms the list for sketch creation.
	 * @param list List of k-mer hash values
	 * @return Sketch created from the list
	 */
	Sketch toSketch(LongList list){
		list.sort();
		assert(list.size==0 || list.get(list.size()-1)>=minHashValue) : list.size+", "+list.get(list.size()-1)+", "+minHashValue;
		list.shrinkToUnique();
		list.reverse();
		for(int i=0; i<list.size; i++){list.array[i]=Long.MAX_VALUE-list.array[i];}
		return new Sketch(list.toArray(), null, null, null, null, null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Recursively traverses k-mer tree and adds qualifying k-mers to heap.
	 * Processes left and right subtrees in depth-first order.
	 * @param heap Target sketch heap
	 * @param kn Current k-mer tree node
	 */
	private void addRecursive(SketchHeap heap, KmerNode kn){
		if(kn==null){return;}
		if(kn.count()>=minKeyOccuranceCount){
			heap.genomeSizeKmers++;
			long kmer=kn.pivot();
			long hashcode=hash(kmer);
			if(hashcode>=minHashValue){heap.add(hashcode);}
		}
		if(kn.left()!=null){addRecursive(heap, kn.left());}
		if(kn.right()!=null){addRecursive(heap, kn.right());}
	}
	
	/**
	 * Recursively traverses k-mer tree and adds qualifying k-mers to list.
	 * Processes left and right subtrees in depth-first order.
	 * @param list Target k-mer list
	 * @param kn Current k-mer tree node
	 */
	private void addRecursive(LongList list, KmerNode kn){
		if(kn==null){return;}
		if(kn.count()>=minKeyOccuranceCount){
			long kmer=kn.pivot();
			long hashcode=hash(kmer);
			if(hashcode>=minHashValue){list.add(hashcode);}
		}
		if(kn.left()!=null){addRecursive(list, kn.left());}
		if(kn.right()!=null){addRecursive(list, kn.right());}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             I/O              ----------------*/
	/*--------------------------------------------------------------*/
	
//	public static ArrayList<Sketch> loadSketches_ST(String...fnames){
//		ArrayList<Sketch> sketches=null;
//		for(String s : fnames){
//			ArrayList<Sketch> temp;
//			if(s.indexOf(',')<0 || s.startsWith("stdin") || new File(s).exists()){
//				temp=loadSketches(s);
//			}else{
//				temp=loadSketches_ST(s.split(","));
//			}
//			if(sketches==null){sketches=temp;}
//			else{sketches.addAll(temp);}
//		}
//		return sketches;
//	}
	
//	public static ArrayList<Sketch> loadSketches_MT(ArrayList<String> fnames){
//		return loadSketches_MT(0, null, fnames.toArray(new String[0]));
//	}
	
	/**
	 * Multithreaded sketch loading from collection of filenames using DisplayParams.
	 * @param params Display parameters for processing configuration
	 * @param fnames Collection of filenames to load sketches from
	 * @return List of loaded sketches
	 */
	public ArrayList<Sketch> loadSketches_MT(DisplayParams params, Collection<String> fnames){
		return loadSketches_MT(params.mode, params.samplerate, params.maxReads,
				params.minEntropy, params.minProb, params.minQual, fnames);
	}
	
	/**
	 * Multithreaded sketch loading from filenames array using DisplayParams.
	 * @param params Display parameters for processing configuration
	 * @param fnames Filenames to load sketches from
	 * @return List of loaded sketches
	 */
	public ArrayList<Sketch> loadSketches_MT(DisplayParams params, String...fnames){
		return loadSketches_MT(params.mode, params.samplerate, params.maxReads,
				params.minEntropy, params.minProb, params.minQual, fnames);
	}
	
	/**
	 * Multithreaded sketch loading from collection with explicit parameters.
	 *
	 * @param mode Processing mode (per file, per sequence, etc.)
	 * @param samplerate Fraction of reads to sample (0.0 to 1.0)
	 * @param reads Maximum number of reads to process
	 * @param minEntropy Minimum sequence entropy threshold
	 * @param minProb Minimum base call probability
	 * @param minQual Minimum quality score threshold
	 * @param fnames Collection of filenames to process
	 * @return List of loaded sketches
	 */
	public ArrayList<Sketch> loadSketches_MT(int mode, float samplerate, long reads, float minEntropy, float minProb, byte minQual, Collection<String> fnames){
		return loadSketches_MT(mode, samplerate, reads, minEntropy, minProb, minQual, fnames.toArray(new String[0]));
	}
	
	//TODO: This is only multithreaded per file in persequence mode.
	/**
	 * Main multithreaded sketch loading method.
	 * Handles comma-separated filename expansion and parallel processing.
	 *
	 * @param mode Processing mode (per file, per sequence, etc.)
	 * @param samplerate Fraction of reads to sample (0.0 to 1.0)
	 * @param reads Maximum number of reads to process
	 * @param minEntropy Minimum sequence entropy threshold
	 * @param minProb Minimum base call probability
	 * @param minQual Minimum quality score threshold
	 * @param fnames Filenames to process (may contain comma-separated lists)
	 * @return List of loaded sketches from all files
	 */
	public ArrayList<Sketch> loadSketches_MT(int mode, float samplerate, long reads, float minEntropy, float minProb, byte minQual, String...fnames){
		
		ConcurrentLinkedQueue<StringNum> decomposedFnames=new ConcurrentLinkedQueue<StringNum>();
		int num=0;
		for(String s : fnames){
			if(s.indexOf(',')<0 || s.startsWith("stdin") || new File(s).exists()){
				num++;
				decomposedFnames.add(new StringNum(s, num));
			}else{
				for(String s2 : s.split(",")){
					num++;
					decomposedFnames.add(new StringNum(s2, num));
				}
			}
		}

		if(decomposedFnames.size()==0){return null;}
		if(decomposedFnames.size()==1){return loadSketchesFromFile(decomposedFnames.poll().s, null, 0, reads, mode, samplerate, minEntropy, minProb, minQual, false);}
		
		//Determine how many threads may be used
		final int threads=Tools.min(Shared.threads(), decomposedFnames.size());
		
		//Fill a list with LoadThreads
		ArrayList<LoadThread> allt=new ArrayList<LoadThread>(threads);
		
		for(int i=0; i<threads; i++){
			allt.add(new LoadThread(decomposedFnames, mode, samplerate, reads, minEntropy, minProb, minQual));
		}
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		
		//Start the threads
		for(LoadThread lt : allt){lt.start();}

		//Wait for completion of all threads
		boolean success=true;
		for(LoadThread lt : allt){

			//Wait until this thread has terminated
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					lt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			sketches.addAll(lt.list);
			success&=lt.success;
		}
		assert(success) : "Failure loading some files.";
		return sketches;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads sketches from a single file using DisplayParams configuration.
	 *
	 * @param fname0 Filename to process
	 * @param smm Optional pre-configured SketchMakerMini instance
	 * @param maxThreads Maximum threads for multithreaded sequence processing
	 * @param reads Maximum number of reads to process
	 * @param mode Processing mode
	 * @param params Display parameters for configuration
	 * @param allowZeroSizeSketch Whether to allow empty sketch creation
	 * @return List of loaded sketches
	 */
	public ArrayList<Sketch> loadSketchesFromFile(final String fname0, SketchMakerMini smm, 
			int maxThreads, long reads, int mode, DisplayParams params, boolean allowZeroSizeSketch){
		return loadSketchesFromFile(fname0, smm, maxThreads, reads, mode,
				params.samplerate, params.minEntropy, params.minProb, params.minQual, allowZeroSizeSketch);
	}
	
	/**
	 * Main file loading method that handles both sequence and sketch files.
	 * Automatically detects file type and routes to appropriate parser.
	 *
	 * @param fname0 Filename to process
	 * @param smm Optional pre-configured SketchMakerMini instance
	 * @param maxThreads Maximum threads for multithreaded sequence processing
	 * @param reads Maximum number of reads to process
	 * @param mode Processing mode (per file, per sequence, etc.)
	 * @param samplerate Fraction of reads to sample
	 * @param minEntropy Minimum sequence entropy threshold
	 * @param minProb Minimum base call probability
	 * @param minQual Minimum quality score threshold
	 * @param allowZeroSizeSketch Whether to allow empty sketch creation
	 * @return List of loaded sketches
	 */
	public ArrayList<Sketch> loadSketchesFromFile(final String fname0, SketchMakerMini smm, 
			int maxThreads, long reads, int mode, float samplerate, float minEntropy, float minProb, byte minQual, boolean allowZeroSizeSketch){
		assert(fname0!=null);//123
		if(fname0==null){return null;}
		FileFormat ff=FileFormat.testInput(fname0, FileFormat.FASTA, null, false, true);
		if(ff.isSequence()){
			return loadSketchesFromSequenceFile(ff, smm, maxThreads, reads, mode, samplerate, minEntropy, minProb, minQual, allowZeroSizeSketch);
		}else{
			return SketchObject.LOADER2 ? loadSketchesFromSketchFile2(ff, allowZeroSizeSketch) : loadSketchesFromSketchFile(ff, allowZeroSizeSketch);
		}
	}
	
	/**
	 * Loads sketches from sequence files (FASTA, FASTQ, SAM).
	 * Uses multithreaded processing for large FASTQ files when beneficial.
	 *
	 * @param ff File format descriptor for the sequence file
	 * @param smm Optional pre-configured SketchMakerMini instance
	 * @param maxThreads Maximum threads for processing
	 * @param reads Maximum number of reads to process
	 * @param mode Processing mode
	 * @param samplerate Fraction of reads to sample
	 * @param minEntropy Minimum sequence entropy threshold
	 * @param minProb Minimum base call probability
	 * @param minQual Minimum quality score threshold
	 * @param allowZeroSizeSketch Whether to allow empty sketch creation
	 * @return List of loaded sketches
	 */
	private ArrayList<Sketch> loadSketchesFromSequenceFile(final FileFormat ff, SketchMakerMini smm, 
			int maxThreads, long reads, int mode, float samplerate, float minEntropy, float minProb, byte minQual, boolean allowZeroSizeSketch){
		maxThreads=(maxThreads<1 ? Shared.threads() : Tools.min(maxThreads, Shared.threads()));
		
//		assert(false) : (ff.fasta() || ff.fastq() || ff.samOrBam())+", "+ff.fastq()+", "+maxThreads+", "+
//				allowMultithreadedFastq+", "+forceDisableMultithreadedFastq+", "+(mode==ONE_SKETCH);
		
		if(ff.fastq() && allowMultithreadedFastq && !forceDisableMultithreadedFastq && (mode==ONE_SKETCH || mode==PER_FILE) &&
				maxThreads>1 && Shared.threads()>2 && (reads<1 || reads*samplerate*(mergePairs ? 2 : 1)>=1000)){//limit is low due to PacBio long reads
			
			final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
			Read.VALIDATE_IN_CONSTRUCTOR=false;

			if(verbose2){System.err.println("Loading a sketch multithreaded.");}
			Sketch sketch=processReadsMT(ff.name(), maxThreads, reads, mode, samplerate, minEntropy, minProb, minQual, allowZeroSizeSketch);
			
			Read.VALIDATE_IN_CONSTRUCTOR=vic;
			
			ArrayList<Sketch> list=new ArrayList<Sketch>(1);
			if(sketch!=null && (sketch.length()>0 || allowZeroSizeSketch)){
				sketch.loadSSU();
				list.add(sketch);
			}
			return list;
		}
		if(smm==null){smm=new SketchMakerMini(this, mode, minEntropy, minProb, minQual);}
		if(verbose2){System.err.println("Loading sketches via SMM.");}
		ArrayList<Sketch> sketches=smm.toSketches(ff.name(), samplerate, reads);
		if(verbose2){System.err.println("Loaded "+(sketches==null ? 0 : sketches.size())+" sketches via SMM.");}
		return sketches;
	}
	
	/**
	 * Loads sketches from text sketch files using line-by-line parsing.
	 * Supports multiple encoding formats (A48, HEX, RAW) and metadata.
	 *
	 * @param ff File format descriptor for the sketch file
	 * @param allowZeroSizeSketch Whether to allow empty sketch creation
	 * @return List of loaded sketches from the file
	 */
	private ArrayList<Sketch> loadSketchesFromSketchFile(final FileFormat ff, boolean allowZeroSizeSketch){
		
		boolean A48=(Sketch.CODING==Sketch.A48), HEX=(Sketch.CODING==Sketch.HEX), NUC=false, delta=true, counts=false, unsorted=false;
		
		if(verbose2){System.err.println("Loading sketches from text.");}
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		ByteFile bf=ByteFile.makeByteFile(ff);
		int currentSketchSize=stTargetSketchSize;
		int taxID=-1;
		long spid=-1;
		long imgID=-1;
		long genomeSizeBases=0, genomeSizeKmers=0, genomeSequences=0;
		long[] baseCounts=null;
		byte[] r16S=null;
		byte[] r18S=null;
		int r16SLen=0;
		int r18SLen=0;
		float probCorrect=-1;
		String name=null, name0=null, fname=ff.simpleName();
		LongList list=null;
		IntList countList=null;
		ArrayList<String> meta=null;
		long sum=0;
		byte[] lastHeader=null;
		
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0){
				//				System.err.println("Processing line "+new String(line));
				if(line[0]=='#'){
					if(r16SLen>0 || r18SLen>0){
						byte[] ssu=KillSwitch.copyOfRange(line, 5, line.length);
						if(Tools.startsWith(line, "#16S:") || Tools.startsWith(line, "#SSU:")){
							assert(r16SLen>0);
							assert(ssu.length==r16SLen) : r16SLen+", "+line.length+"\n"+new String(line)+"\n";
							r16S=ssu;
							r16SLen=0;
						}else if(Tools.startsWith(line, "#18S:")){
							assert(r18SLen>0);
							assert(ssu.length==r18SLen) : r18SLen+", "+line.length+"\n"+new String(line)+"\n";
							r18S=ssu;
							r18SLen=0;
						}else{
							assert(false) : new String(line);
						}
					}else{
						lastHeader=line;
						if(list!=null){
							assert(list.size==list.array.length);
							if(NUC || unsorted){
								list.sort();
								list.shrinkToUnique();
							}else{
								list.shrink();
							}
							if(list.size()>0 || allowZeroSizeSketch){
								int[] keyCounts=countList==null ? null : countList.array;
								Sketch sketch=new Sketch(list.array, keyCounts, baseCounts, r16S, r18S, taxID, imgID, 
										genomeSizeBases, genomeSizeKmers, genomeSequences, probCorrect, name, name0, fname, meta);
								sketch.spid=spid;
								sketches.add(sketch);
							}
							//						System.err.println("Made sketch "+sketch);
						}
						name=name0=null;
						fname=ff.simpleName();
						list=null;
						countList=null;
						meta=null;
						baseCounts=null;
						r16S=null;
						r18S=null;
						r16SLen=0;
						r18SLen=0;
						sum=0;
						taxID=-1;
						imgID=-1;
						genomeSizeBases=0;
						genomeSizeKmers=0;
						genomeSequences=0;
						probCorrect=-1;
						int k_sketch=defaultK;
						int k2_sketch=0;
						int hashVersion_sketch=1;

						if(line.length>1){
							String[] split=new String(line, 1, line.length-1).split("\t");
							for(String s : split){
								final int colon=s.indexOf(':');
								final String sub=s.substring(colon+1);
								if(s.startsWith("SZ:") || s.startsWith("SIZE:")){//Sketch length
									currentSketchSize=Integer.parseInt(sub);
								}else if(s.startsWith("CD:")){//Coding
									A48=HEX=NUC=delta=counts=unsorted=false;

									for(int i=0; i<sub.length(); i++){
										char c=sub.charAt(i);
										if(c=='A'){A48=true;}
										else if(c=='H'){HEX=true;}
										else if(c=='R'){A48=HEX=false;}
										else if(c=='N'){NUC=true;}
										else if(c=='D'){delta=true;}
										else if(c=='C'){counts=true;}
										else if(c=='U'){unsorted=true;}
										else if(c=='M'){assert(aminoOrTranslate()) : "Amino sketch in non-amino mode: "+new String(line);}
										else if(c=='8'){assert(amino8) : "Amino8 sketch in non-amino8 mode: "+new String(line);}
										else{assert(false) : "Unknown coding symbol: "+c+"\t"+new String(line);}
									}
								}else if(s.startsWith("K:")){//Kmer length
									if(sub.indexOf(',')>=0){
										String[] subsplit=sub.split(",");
										assert(subsplit.length==2) : "Bad header component "+s;
										int x=Integer.parseInt(subsplit[0]);
										int y=Integer.parseInt(subsplit[1]);
										k_sketch=Tools.max(x, y);
										k2_sketch=Tools.min(x, y);
									}else{
										k_sketch=Integer.parseInt(sub);
										k2_sketch=0;
									}
								}else if(s.startsWith("H:")){//Hash version
									hashVersion_sketch=Integer.parseInt(sub);
								}else if(s.startsWith("BC:") || s.startsWith("BASECOUNTS:")){//ACGTN
									baseCounts=Parse.parseLongArray(sub);
								}else if(s.startsWith("GS:") || s.startsWith("GSIZE:")){//Genomic bases
									genomeSizeBases=Long.parseLong(sub);
								}else if(s.startsWith("GK:") || s.startsWith("GKMERS:")){//Genomic kmers
									genomeSizeKmers=Long.parseLong(sub);
								}else if(s.startsWith("GQ:")){
									genomeSequences=Long.parseLong(sub);
								}else if(s.startsWith("GE:")){//Genome size estimate kmers
									//ignore
								}else if(s.startsWith("PC:")){//Probability of correctness
									probCorrect=Float.parseFloat(sub);
								}else if(s.startsWith("ID:") || s.startsWith("TAXID:")){
									taxID=Integer.parseInt(sub);
								}else if(s.startsWith("IMG:")){
									imgID=Long.parseLong(sub);
								}else if(s.startsWith("SPID:")){
									spid=Integer.parseInt(sub);
								}else if(s.startsWith("NM:") || s.startsWith("NAME:")){
									name=sub;
								}else if(s.startsWith("FN:")){
									fname=sub;
								}else if(s.startsWith("NM0:")){
									name0=sub;
								}else if(s.startsWith("MT_")){
									if(meta==null){meta=new ArrayList<String>(1);}
									meta.add(s.substring(3));
								}else if(s.startsWith("16S:") || s.startsWith("SSU:")){
									r16SLen=Integer.parseInt(sub);
								}else if(s.startsWith("18S:")){
									r18SLen=Integer.parseInt(sub);
								}else{
									assert(false) : "Unsupported header tag "+s;
								}
							}
						}

						if(KILL_OK){
							if(k_sketch!=k && !NUC){KillSwitch.kill("Sketch kmer length "+k_sketch+
									" differs from loaded kmer length "+k+"\n"+new String(line)+"\nfile: "+ff.name());}
							if(k2_sketch!=k2 && !NUC){KillSwitch.kill("Sketch kmer length "+k_sketch+","+k2_sketch+
									" differs from loaded kmer length "+k+","+k2+"\n"+new String(line)+"\nfile: "+ff.name());}
							if(hashVersion_sketch!=HASH_VERSION && !NUC){KillSwitch.kill("Sketch hash version "+hashVersion_sketch+
									" differs from loaded hash version "+HASH_VERSION+".\n"
									+ "You may need to download the latest version of BBTools.\n"+new String(line)+"\nfile: "+ff.name());}
						}else{//Potential hang
							assert(k_sketch==k || NUC) : "Sketch kmer length "+k_sketch+" differs from loaded kmer length "+k+"\n"+new String(line);
							assert(k2_sketch==k2 || NUC) : "Sketch kmer length "+k_sketch+","+k2_sketch+" differs from loaded kmer length "+k+","+k2+"\n"+new String(line);
							assert(hashVersion_sketch==HASH_VERSION || NUC) : "Sketch hash version "+hashVersion_sketch+
							" differs from loaded hash version "+HASH_VERSION+".\n"
							+ "You may need to download the latest version of BBTools.\n"+new String(line)+"\n";
						}
						if(currentSketchSize>0 || allowZeroSizeSketch){
							list=new LongList(Tools.max(1, currentSketchSize));
							if(counts){countList=new IntList(Tools.max(1, currentSketchSize));}
						}
					}
				}else{
					long x=(counts ? Sketch.parseA48C(line, countList) : A48 ? Sketch.parseA48(line) :
						HEX ? Sketch.parseHex(line) : NUC ? Sketch.parseNuc(line) : Parse.parseLong(line));
//					System.err.println("sum="+sum+", x="+x+" -> "+(sum+x));
					sum+=x;
					assert(x>=0 || NUC) : "x="+x+"\nline="+new String(line)+"\nheader="+(lastHeader==null ? "null" : new String(lastHeader))+"\nlineNum="+bf.lineNum()+"\n";
					assert(sum>=0 || !delta) : "The sketch was made with delta compression off.  Please regenerate it.";
					assert(list!=null) : new String(line);
					long key=(delta ? sum : x);
					if(key>=0){list.add(key);}
				}
			}
		}
		if(list!=null && (list.size>0 || allowZeroSizeSketch)){
			assert(list.size==list.array.length || NUC || unsorted || (allowZeroSizeSketch && list.size==0)) : list.size+"!="+list.array.length;
			if(NUC || unsorted){
				list.sort();
				list.shrinkToUnique();
			}else{
				list.shrink();
			}
			int[] keyCounts=countList==null ? null : countList.array;
			Sketch sketch=new Sketch(list.array, keyCounts, baseCounts, r16S, r18S, taxID, imgID, 
					genomeSizeBases, genomeSizeKmers, genomeSequences, probCorrect, name, name0, fname, meta);
			sketch.spid=spid;
			sketches.add(sketch);
		}
		if(verbose2){System.err.println("Loaded "+sketches.size()+" sketches from text.");}
		return sketches;
	}
	
	/** Usually much faster due to not manifesting the multithreaded load Java slowdown.  Should incur less garbage collection also. */
	private ArrayList<Sketch> loadSketchesFromSketchFile2(final FileFormat ff, boolean allowZeroSizeSketch){
		
		boolean A48=(Sketch.CODING==Sketch.A48), HEX=(Sketch.CODING==Sketch.HEX), NUC=false, delta=true, counts=false, unsorted=false;
		
		if(verbose2){System.err.println("Loading sketches from text.");}
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		
		InputStream is=ReadWrite.getInputStream(ff.name(), BUFFERED_READER, false, false);
		byte[] buffer=new byte[BUFLEN];
		int start, limit=0;
		try {
			limit=is.read(buffer);
		} catch (IOException e) {
			KillSwitch.exceptionKill(e);
		}
		
		
		int currentSketchSize=stTargetSketchSize;
		int taxID=-1;
		long spid=-1;
		long imgID=-1;
		long genomeSizeBases=0, genomeSizeKmers=0, genomeSequences=0;
		long[] baseCounts=null;
		byte[] r16S=null;
		byte[] r18S=null;
		int r16SLen=0;
		int r18SLen=0;
		float probCorrect=-1;
		String name=null, name0=null, fname=ff.simpleName();
		LongList list=null;
		IntList countList=null;
		ArrayList<String> meta=null;
		long sum=0;
		byte[] lastHeader=null;
		ByteBuilder bb=new ByteBuilder(256);
		
		for(start=0; start<limit;){
			//				System.err.println("Processing line "+new String(line));
			if(buffer[start]=='#'){
				bb.clear();
				try {
					while(buffer[start]!='\n'){
						bb.append(buffer[start]);
						start++;
						if(start>=limit){start=0; limit=is.read(buffer);}
					}
				} catch (IOException e) {
					KillSwitch.exceptionKill(e);
				}
				start++;
				
				if(r16SLen>0 || r18SLen>0){
					byte[] ssu=bb.toBytes(5, bb.length());
					if(bb.startsWith("#16S:") || bb.startsWith("#SSU:")){
						assert(r16SLen>0);
						assert(ssu.length==r16SLen) : r16SLen+", "+bb.length+"\n"+bb+"\n";
						r16S=ssu;
						r16SLen=0;
					}else if(bb.startsWith("#18S:")){
						assert(r18SLen>0);
						assert(ssu.length==r18SLen) : r18SLen+", "+bb.length+"\n"+bb+"\n";
						r18S=ssu;
						r18SLen=0;
					}else{
						assert(false) : bb;
					}
				}else{

					//				byte[] line=lastHeader=bb.toBytes();
					if(list!=null){

						//This assertion fails sometimes for Silva per-sequence mode, but it's not important
						//					assert(list.size==list.array.length) : list.size+", "+list.array.length+(lastHeader==null ? "" : ", "+new String(lastHeader));

						if(NUC || unsorted){
							list.sort();
							list.shrinkToUnique();
						}else{
							list.shrink();
						}
						if(list.size()>0 || allowZeroSizeSketch){
							int[] keyCounts=countList==null ? null : countList.array;
							Sketch sketch=new Sketch(list.array, keyCounts, baseCounts, r16S, r18S, taxID, imgID, 
									genomeSizeBases, genomeSizeKmers, genomeSequences, probCorrect, name, name0, fname, meta);
							sketch.spid=spid;
							sketches.add(sketch);
						}
						//						System.err.println("Made sketch "+sketch);
					}
					name=name0=null;
					fname=ff.simpleName();
					list=null;
					countList=null;
					meta=null;
					baseCounts=null;
					r16S=null;
					r18S=null;
					r16SLen=0;
					r18SLen=0;
					sum=0;
					taxID=-1;
					imgID=-1;
					genomeSizeBases=0;
					genomeSizeKmers=0;
					genomeSequences=0;
					probCorrect=-1;
					int k_sketch=defaultK;
					int k2_sketch=0;
					int hashVersion_sketch=1;

					if(bb.length>1){
						String[] split=new String(bb.array, 1, bb.length-1).split("\t");
						for(String s : split){
							final int colon=s.indexOf(':');
							final String sub=s.substring(colon+1);
							if(s.startsWith("SZ:") || s.startsWith("SIZE:")){//Sketch length
								currentSketchSize=Integer.parseInt(sub);
							}else if(s.startsWith("CD:")){//Coding
								A48=HEX=NUC=delta=counts=unsorted=false;

								for(int i=0; i<sub.length(); i++){
									char c=sub.charAt(i);
									if(c=='A'){A48=true;}
									else if(c=='H'){HEX=true;}
									else if(c=='R'){A48=HEX=false;}
									else if(c=='N'){NUC=true;}
									else if(c=='D'){delta=true;}
									else if(c=='C'){counts=true;}
									else if(c=='U'){unsorted=true;}
									else if(c=='M'){assert(aminoOrTranslate()) : "Amino sketch in non-amino mode: "+bb;}
									else if(c=='8'){assert(amino8) : "Amino8 sketch in non-amino8 mode: "+bb;}
									else{assert(false) : "Unknown coding symbol: "+c+"\t"+bb;}
								}
							}else if(s.startsWith("K:")){//Kmer length
								if(sub.indexOf(',')>=0){
									String[] subsplit=sub.split(",");
									assert(subsplit.length==2) : "Bad header component "+s;
									int x=Integer.parseInt(subsplit[0]);
									int y=Integer.parseInt(subsplit[1]);
									k_sketch=Tools.max(x, y);
									k2_sketch=Tools.min(x, y);
								}else{
									k_sketch=Integer.parseInt(sub);
									k2_sketch=0;
								}
							}else if(s.startsWith("H:")){//Hash version
								hashVersion_sketch=Integer.parseInt(sub);
							}else if(s.startsWith("BC:") || s.startsWith("BASECOUNTS:")){//ACGTN
								baseCounts=Parse.parseLongArray(sub);
							}else if(s.startsWith("GS:") || s.startsWith("GSIZE:")){//Genomic bases
								genomeSizeBases=Long.parseLong(sub);
							}else if(s.startsWith("GK:") || s.startsWith("GKMERS:")){//Genomic kmers
								genomeSizeKmers=Long.parseLong(sub);
							}else if(s.startsWith("GQ:")){
								genomeSequences=Long.parseLong(sub);
							}else if(s.startsWith("GE:")){//Genome size estimate kmers
								//ignore
							}else if(s.startsWith("PC:")){//Probability of correctness
								probCorrect=Float.parseFloat(sub);
							}else if(s.startsWith("ID:") || s.startsWith("TAXID:")){
								taxID=Integer.parseInt(sub);
							}else if(s.startsWith("IMG:")){
								imgID=Long.parseLong(sub);
							}else if(s.startsWith("SPID:")){
								spid=Integer.parseInt(sub);
							}else if(s.startsWith("NM:") || s.startsWith("NAME:")){
								name=sub;
							}else if(s.startsWith("FN:")){
								fname=sub;
							}else if(s.startsWith("NM0:")){
								name0=sub;
							}else if(s.startsWith("MT_")){
								if(meta==null){meta=new ArrayList<String>(1);}
								meta.add(s.substring(3));
							}else if(s.startsWith("16S:") || s.startsWith("SSU:")){
								r16SLen=Integer.parseInt(sub);
							}else if(s.startsWith("18S:")){
								r18SLen=Integer.parseInt(sub);
							}else{
								assert(false) : "Unsupported header tag "+s;
							}
						}
					}

					if(KILL_OK){
						if(k_sketch!=k && !NUC){KillSwitch.kill("Sketch kmer length "+k_sketch+
								" differs from loaded kmer length "+k+"\n"+bb+"\nfile: "+ff.name());}
						if(k2_sketch!=k2 && !NUC){KillSwitch.kill("Sketch kmer length "+k_sketch+","+k2_sketch+
								" differs from loaded kmer length "+k+","+k2+"\n"+bb+"\nfile: "+ff.name());}
						if(hashVersion_sketch!=HASH_VERSION && !NUC){KillSwitch.kill("Sketch hash version "+hashVersion_sketch+
								" differs from loaded hash version "+HASH_VERSION+".\n"
								+ "You may need to download the latest version of BBTools.\nfile: "+ff.name());}
					}else{//Potential hang
						assert(k_sketch==k || NUC) : "Sketch kmer length "+k_sketch+" differs from loaded kmer length "+k+"\n"+bb;
						assert(k2_sketch==k2 || NUC) : "Sketch kmer length "+k_sketch+","+k2_sketch+" differs from loaded kmer length "+k+","+k2+"\n"+bb;
						assert(hashVersion_sketch==HASH_VERSION || NUC) : "Sketch hash version "+hashVersion_sketch+
						" differs from loaded hash version "+HASH_VERSION+".\n"
						+ "You may need to download the latest version of BBTools.\nfile: "+ff.name();
					}
					if(currentSketchSize>0 || allowZeroSizeSketch){
						list=new LongList(Tools.max(1, currentSketchSize));
						if(counts){countList=new IntList(Tools.max(1, currentSketchSize));}
					}
				}
			}else{
				bb.clear();
				try {
					while(buffer[start]!='\n'){
						bb.append(buffer[start]);
						start++;
						if(start>=limit){start=0; limit=is.read(buffer);}
					}
				} catch (IOException e) {
					KillSwitch.exceptionKill(e);
				}
				start++;
				
				long x=(counts ? Sketch.parseA48C(bb, countList) : A48 ? Sketch.parseA48(bb) :
					HEX ? Sketch.parseHex(bb) : NUC ? Sketch.parseNuc(bb) : Parse.parseLong(bb.array, 0, bb.length));
				//					System.err.println("sum="+sum+", x="+x+" -> "+(sum+x));
				sum+=x;
				assert(x>=0 || NUC) : "x="+x+"\nline="+bb+"\nheader="+(lastHeader==null ? "null" : new String(lastHeader))+"\n";
				assert(sum>=0 || !delta) : "The sketch was made with delta compression off.  Please regenerate it.";
				assert(list!=null) : bb;
				long key=(delta ? sum : x);
				if(key>=0){list.add(key);}
			}
			
			if(start>=limit){
				start=0;
				try {
					limit=is.read(buffer);
				} catch (IOException e) {
					KillSwitch.exceptionKill(e);
				}
			}
		}
		if(list!=null && (list.size>0 || allowZeroSizeSketch)){
			
			//This assertion fails sometimes for Silva per-sequence mode, but it's not important
//			assert(list.size==list.array.length || NUC || unsorted || (allowZeroSizeSketch && list.size==0)) : 
//				list.size+"!="+list.array.length+(lastHeader==null ? "" : "\n"+new String(lastHeader));
			
			if(NUC || unsorted){
				list.sort();
				list.shrinkToUnique();
			}else{
				list.shrink();
			}
			int[] keyCounts=countList==null ? null : countList.array;
			Sketch sketch=new Sketch(list.array, keyCounts, baseCounts, r16S, r18S, taxID, imgID, 
					genomeSizeBases, genomeSizeKmers, genomeSequences, probCorrect, name, name0, fname, meta);
			sketch.spid=spid;
			sketches.add(sketch);
		}
		if(verbose2){System.err.println("Loaded "+sketches.size()+" sketches from text.");}
		try {
			is.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sketches;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads sketches from a string containing sketch data.
	 * Parses the same format as sketch files but from memory.
	 * @param sketchString String containing sketch data in text format
	 * @return List of loaded sketches
	 */
	public ArrayList<Sketch> loadSketchesFromString(String sketchString){
		boolean A48=(Sketch.CODING==Sketch.A48), HEX=(Sketch.CODING==Sketch.HEX), NUC=false, delta=true, counts=false, unsorted=false;
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		int currentSketchSize=stTargetSketchSize;
		int taxID=-1;
		long spid=-1;
		long imgID=-1;
		long genomeSizeBases=0, genomeSizeKmers=0, genomeSequences=0;
		long[] baseCounts=null;
		byte[] r16S=null;
		byte[] r18S=null;
		int r16SLen=0;
		int r18SLen=0;
		float probCorrect=-1;
		String name=null, name0=null, fname=null;
		LongList list=null;
		IntList countList=null;
		ArrayList<String> meta=null;
		long sum=0;
		
		String[] split0=sketchString.split("\n");
		for(String line : split0){
			if(line.length()>0){
//				System.err.println("Processing line "+new String(line));
				if(line.charAt(0)=='#'){
					if(line.length()>1 && line.charAt(1)=='#'){
						//ignore
					}else if(r16SLen>0 || r18SLen>0){
						byte[] ssu=KillSwitch.copyOfRange(line.getBytes(), 5, line.length());
						if(line.startsWith("#16S:") || line.startsWith("#SSU:")){
							assert(r16SLen>0);
							assert(ssu.length==r16SLen) : r16SLen+", "+line.length()+"\n"+line+"\n";
							r16S=ssu;
							r16SLen=0;
						}else if(line.startsWith("#18S:")){
							assert(r18SLen>0);
							assert(ssu.length==r18SLen) : r18SLen+", "+line.length()+"\n"+line+"\n";
							r18S=ssu;
							r18SLen=0;
						}else{
							assert(false) : line;
						}
					}else{
						if(list!=null){
							assert(list.size==list.array.length);
							if(NUC || unsorted){
								list.sort();
								list.shrinkToUnique();
							}else{
								list.shrink();
							}
							if(list.size()>0){
								int[] keyCounts=countList==null ? null : countList.array;
								Sketch sketch=new Sketch(list.array, keyCounts, baseCounts, r16S, r18S, taxID, imgID, 
										genomeSizeBases, genomeSizeKmers, genomeSequences, probCorrect, name, name0, fname, meta);
								sketch.spid=spid;
								sketches.add(sketch);
							}
							//						System.err.println("Made sketch "+sketch);
						}
						name=name0=null;
						fname=null;
						list=null;
						countList=null;
						meta=null;
						baseCounts=null;
						r16S=null;
						r18S=null;
						r16SLen=0;
						r18SLen=0;
						sum=0;
						taxID=-1;
						imgID=-1;
						genomeSizeBases=0;
						genomeSizeKmers=0;
						genomeSequences=0;
						probCorrect=-1;
						int k_sketch=defaultK;
						int k2_sketch=0;
						int hashVersion_sketch=1;

						if(line.length()>1){
							String[] split=line.substring(1).split("\t");
							for(String s : split){
								final int colon=s.indexOf(':');
								final String sub=s.substring(colon+1);
								if(s.startsWith("SZ:") || s.startsWith("SIZE:")){//Sketch length
									currentSketchSize=Integer.parseInt(sub);
								}else if(s.startsWith("CD:")){//Coding
									A48=HEX=NUC=delta=counts=false;
									
									for(int i=0; i<sub.length(); i++){
										char c=sub.charAt(i);
										if(c=='A'){A48=true;}
										else if(c=='H'){HEX=true;}
										else if(c=='R'){A48=HEX=false;}
										else if(c=='N'){NUC=true;}
										else if(c=='D'){delta=true;}
										else if(c=='C'){counts=true;}
										else if(c=='U'){unsorted=true;}
										else if(c=='M'){assert(aminoOrTranslate()) : "Amino sketch in non-amino mode: "+new String(line);}
										else if(c=='8'){assert(amino8) : "Amino8 sketch in non-amino8 mode: "+new String(line);}
										else{assert(false) : "Unknown coding symbol: "+c+"\t"+new String(line);}
									}
									
								}else if(s.startsWith("K:")){//Kmer length
									if(sub.indexOf(',')>=0){
										String[] subsplit=sub.split(",");
										assert(subsplit.length==2) : "Bad header component "+s;
										int x=Integer.parseInt(subsplit[0]);
										int y=Integer.parseInt(subsplit[1]);
										k_sketch=Tools.max(x, y);
										k2_sketch=Tools.min(x, y);
									}else{
										k_sketch=Integer.parseInt(s);
										k2_sketch=0;
									}
								}else if(s.startsWith("H:")){//Hash version
									hashVersion_sketch=Integer.parseInt(sub);
								}else if(s.startsWith("BC:") || s.startsWith("BASECOUNTS:")){//ACGTN
									baseCounts=Parse.parseLongArray(sub);
								}else if(s.startsWith("GS:") || s.startsWith("GSIZE:")){//Genomic bases
									genomeSizeBases=Long.parseLong(sub);
								}else if(s.startsWith("GK:") || s.startsWith("GKMERS:")){//Genomic kmers
									genomeSizeKmers=Long.parseLong(sub);
								}else if(s.startsWith("GQ:")){
									genomeSequences=Long.parseLong(sub);
								}else if(s.startsWith("GE:")){//Genome size estimate kmers
									//ignore
								}else if(s.startsWith("PC:")){//Probability of correctness
									probCorrect=Float.parseFloat(sub);
								}else if(s.startsWith("ID:") || s.startsWith("TAXID:")){
									taxID=Integer.parseInt(sub);
								}else if(s.startsWith("IMG:")){
									imgID=Long.parseLong(sub);
								}else if(s.startsWith("SPID:")){
									spid=Integer.parseInt(sub);
								}else if(s.startsWith("NM:") || s.startsWith("NAME:")){
									name=sub;
								}else if(s.startsWith("FN:")){
									fname=sub;
								}else if(s.startsWith("NM0:")){
									name0=sub;
								}else if(s.startsWith("MT_")){
									if(meta==null){meta=new ArrayList<String>(1);}
									meta.add(s.substring(3));
								}else if(s.startsWith("16S:") || s.startsWith("SSU:")){
									r16SLen=Integer.parseInt(sub);
								}else if(s.startsWith("18S:")){
									r18SLen=Integer.parseInt(sub);
								}else{
									assert(false) : "Unsupported header tag "+s;
								}
							}
						}
						if(KILL_OK){
							if(k_sketch!=k && !NUC){KillSwitch.kill("Sketch kmer length "+k_sketch+" differs from loaded kmer length "+k+"\n"+new String(line));}
							if(k2_sketch!=k2 && !NUC){KillSwitch.kill("Sketch kmer length "+k_sketch+","+k2_sketch+" differs from loaded kmer length "+k+","+k2+"\n"+new String(line));}
							if(hashVersion_sketch!=HASH_VERSION && !NUC){KillSwitch.kill("Sketch hash version "+hashVersion_sketch+
									" differs from loaded hash version "+HASH_VERSION+".\n"
											+ "You may need to download the latest version of BBTools.\n"+new String(line)+"\n");}
						}else{//Potential hang
							assert(k_sketch==k && !NUC) : "Sketch kmer length "+k_sketch+" differs from loaded kmer length "+k+"\n"+new String(line);
							assert(k2_sketch==k2 && !NUC) : "Sketch kmer length "+k_sketch+","+k2_sketch+" differs from loaded kmer length "+k+","+k2+"\n"+new String(line);
							assert(hashVersion_sketch==HASH_VERSION || NUC) : "Sketch hash version "+hashVersion_sketch+
									" differs from loaded hash version "+HASH_VERSION+".\n"
											+ "You may need to download the latest version of BBTools.\n"+new String(line)+"\n";
						}
						
						
						if(currentSketchSize>0){
							list=new LongList(currentSketchSize);
							if(counts){countList=new IntList(currentSketchSize);}
						}
					}
				}else{
					long x=(counts ? Sketch.parseA48C(line, countList) : A48 ? Sketch.parseA48(line) :
						HEX ? Sketch.parseHex(line) : NUC ? Sketch.parseNuc(line) : Long.parseLong(line));
//					System.err.println("sum="+sum+", x="+x+" -> "+(sum+x));
					sum+=x;
					assert(x>=0 || NUC) : x+"\n"+new String(line);
					assert(sum>=0 || !delta) : "The sketch was made with delta compression off.  Please regenerate it.";
					long key=(delta ? sum : x);
					if(key>=0){list.add(key);}
				}
			}
		}
		
		if(list!=null){
			assert(list.size==list.array.length || list.size()==0 || NUC || unsorted);
			if(NUC || unsorted){
				list.sort();
				list.shrinkToUnique();
			}else{
				list.shrink();
			}
			int[] keyCounts=countList==null ? null : countList.array;
			Sketch sketch=new Sketch(list.array, keyCounts, baseCounts, r16S, r18S, taxID, imgID, 
					genomeSizeBases, genomeSizeKmers, genomeSequences, probCorrect, name, name0, fname, meta);
			sketch.spid=spid;
			sketches.add(sketch);
		}
		return sketches;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Threads            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	
	/** Worker thread for multithreaded sketch file loading.
	 * Processes files from a shared queue and collects resulting sketches. */
	private class LoadThread extends Thread{
		
		/**
		 * Constructs a LoadThread with processing parameters.
		 *
		 * @param queue_ Shared queue of files to process
		 * @param mode_ Processing mode
		 * @param samplerate_ Sampling rate for reads
		 * @param reads_ Maximum reads to process per file
		 * @param minEntropy Minimum sequence entropy threshold
		 * @param minProb Minimum base call probability
		 * @param minQual Minimum quality score threshold
		 */
		public LoadThread(ConcurrentLinkedQueue<StringNum> queue_, int mode_, float samplerate_, long reads_, float minEntropy, float minProb, byte minQual) {
			queue=queue_;
			list=new ArrayList<Sketch>();
			smm=new SketchMakerMini(SketchTool.this, mode_, minEntropy, minProb, minQual);
			samplerate=samplerate_;
			reads=reads_;
		}
		
		@Override
		public void run(){
			success=false;
			for(StringNum sn=queue.poll(); sn!=null; sn=queue.poll()){
				ArrayList<Sketch> temp=null;
				try {
					temp=loadSketchesFromFile(sn.s, smm, 1, reads, smm.mode, samplerate, smm.minEntropy(), smm.minProb(), smm.minQual(), false);
				} catch (Throwable e) {
					System.err.println("Failure loading "+sn+":\n"+e);
					e.printStackTrace();
					success=false;
				}
				if(temp!=null && temp.size()>0){
					if(smm.mode==PER_FILE){
//						assert(temp.size()==1) : temp.size();
						temp.get(0).sketchID=(int)sn.n;
					}
					for(Sketch s : temp){add(s);}
				}
			}
			success=true;
		}
		
		private void add(Sketch s){
			if(list!=null){
				list.add(s);
				return;
			}
			assert(false) : "Unsupported."; //The map logic is broken; needs to be synchronized.
//			if(s.taxID<0){return;}
////			assert(s.taxID>-1) : s.toHeader();
//			TaxNode tn=tree.getNode(s.taxID);
//			while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
//				TaxNode temp=tree.getNode(tn.pid);
//				if(temp==null){break;}
//				tn=temp;
//			}
//			if(tn==null){return;}
//			Integer key=tn.id;
//			Sketch old=map.get(key);
//			if(old==null){
//				s.taxID=key;
//				map.put(key, s);
//			}else{
//				synchronized(old){
//					old.add(s, maxLen);
//				}
//			}
		}
		
		/** Shared queue of filenames for LoadThread processing */
		final ConcurrentLinkedQueue<StringNum> queue;
		/** List to collect sketches loaded by this thread */
		ArrayList<Sketch> list;
		/** Success flag indicating if thread completed without errors */
		boolean success=false;
		/** SketchMakerMini instance for processing sequences */
		final SketchMakerMini smm;
		/** Sampling rate for read processing */
		final float samplerate;
		/** Maximum number of reads to process */
		final long reads;
		
//		ConcurrentHashMap<Integer, Sketch> map;
		
	}
	
	/** Worker thread for multithreaded sequence file processing.
	 * Processes read lists from a concurrent input stream to generate sketches. */
	private class LoadThread2 extends Thread{
		
		/**
		 * Constructs a LoadThread2 for sequence processing.
		 *
		 * @param cris_ Concurrent read input stream
		 * @param minEntropy Minimum sequence entropy threshold
		 * @param minProb Minimum base call probability
		 * @param minQual Minimum quality score threshold
		 */
		LoadThread2(Streamer cris_, float minEntropy, float minProb, byte minQual){
			cris=cris_;
			smm=new SketchMakerMini(SketchTool.this, ONE_SKETCH, minEntropy, minProb, minQual);
		}
		
		@Override
		public void run(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null){//ln!=null prevents a compiler potential null access warning

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					if(r1!=null && !r1.validated()){r1.validate(true);}
					if(r2!=null && !r2.validated()){r2.validate(true);}
					
					smm.processReadPair(r1, r2);
				}

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}
		
		/** Concurrent read input stream for LoadThread2 */
		Streamer cris;
		/** SketchMakerMini instance for sequence processing */
		SketchMakerMini smm;
		
	}

	/** Converts KmerTableSets to Heaps */
	private class SketchThread extends Thread {

		/**
		 * Constructs a SketchThread for parallel k-mer table processing.
		 * @param next_ Atomic counter for distributing work
		 * @param kts_ K-mer table set to process
		 */
		SketchThread(AtomicInteger next_, KmerTableSet kts_){
			next=next_;
			kts=kts_;
		}

		@Override
		public void run(){
			final int ways=kts.ways();
			int tnum=next.getAndIncrement();
			while(tnum<ways){
				HashArray1D table=kts.getTable(tnum);
				if(stTargetSketchSize>0){
					if(heap==null){heap=new SketchHeap(stTargetSketchSize, minKeyOccuranceCount, trackCounts);}
					toHeap(table, heap);
				}else{
					if(list==null){list=new LongList();}
					toList(table, list);
				}
				tnum=next.getAndIncrement();
			}
		}

		/** Atomic counter for distributing work among SketchThreads */
		final AtomicInteger next;
		/** K-mer table set being processed by this thread */
		final KmerTableSet kts;
		/** Sketch heap for collecting k-mers (when using size limit) */
		SketchHeap heap;
		/** List for collecting k-mers (when no size limit) */
		LongList list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Read Loading         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Multithreaded sketch generation from sequence file reads.
	 * Handles paired-end filename patterns with # placeholder.
	 *
	 * @param fname Input filename (may contain # for paired reads)
	 * @param maxThreads Maximum processing threads
	 * @param reads Maximum reads to process
	 * @param mode Processing mode
	 * @param samplerate Sampling rate for reads
	 * @param minEntropy Minimum entropy threshold
	 * @param minProb Minimum probability threshold
	 * @param minQual Minimum quality threshold
	 * @param allowZeroSizeSketch Whether to allow empty sketches
	 * @return Generated sketch from the reads
	 */
	public Sketch processReadsMT(String fname, int maxThreads, long reads, int mode, float samplerate,
			float minEntropy, float minProb, byte minQual, boolean allowZeroSizeSketch){
		if(fname.indexOf('#')>=0 && FileFormat.isFastqExt(ReadWrite.rawExtension(fname)) && !new File(fname).exists()){
			return processReadsMT(fname.replaceFirst("#", "1"), fname.replaceFirst("#", "2"), maxThreads, reads, mode, samplerate,
					minEntropy, minProb, minQual, allowZeroSizeSketch);
		}else{
			return processReadsMT(fname, null, maxThreads, reads, mode, samplerate, minEntropy, minProb, minQual, allowZeroSizeSketch);
		}
	}
	
	/**
	 * Multithreaded sketch generation from paired sequence files.
	 *
	 * @param fname1 First read file
	 * @param fname2 Second read file (may be null)
	 * @param maxThreads Maximum processing threads
	 * @param reads Maximum reads to process
	 * @param mode Processing mode
	 * @param samplerate Sampling rate for reads
	 * @param minEntropy Minimum entropy threshold
	 * @param minProb Minimum probability threshold
	 * @param minQual Minimum quality threshold
	 * @param allowZeroSizeSketch Whether to allow empty sketches
	 * @return Generated sketch from the reads
	 */
	public Sketch processReadsMT(String fname1, String fname2, int maxThreads, long reads, int mode, float samplerate,
			float minEntropy, float minProb, byte minQual, boolean allowZeroSizeSketch){
		final FileFormat ffin1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
		final FileFormat ffin2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
		return processReadsMT(ffin1, ffin2, maxThreads, reads, mode, samplerate, minEntropy, minProb, minQual, allowZeroSizeSketch);
	}
	
	/**
	 * Multithreaded read processing with DisplayParams configuration.
	 *
	 * @param ffin1 First input file format
	 * @param ffin2 Second input file format (may be null)
	 * @param maxThreads Maximum processing threads
	 * @param reads Maximum reads to process
	 * @param mode Processing mode
	 * @param params Display parameters for configuration
	 * @param allowZeroSizeSketch Whether to allow empty sketches
	 * @return Generated sketch from the reads
	 */
	public Sketch processReadsMT(FileFormat ffin1, FileFormat ffin2, int maxThreads, long reads, int mode, DisplayParams params, boolean allowZeroSizeSketch){
		return processReadsMT(ffin1, ffin2, maxThreads, reads, mode, params.samplerate, params.minEntropy, params.minProb, params.minQual, allowZeroSizeSketch);
	}
	
	/**
	 * Main multithreaded read processing implementation.
	 * Creates concurrent read stream and worker threads for parallel k-mer
	 * extraction and sketch generation.
	 *
	 * @param ffin1 First input file format
	 * @param ffin2 Second input file format (may be null)
	 * @param maxThreads Maximum processing threads
	 * @param reads Maximum reads to process
	 * @param mode Processing mode (must be ONE_SKETCH or PER_FILE)
	 * @param samplerate Sampling rate for reads
	 * @param minEntropy Minimum entropy threshold
	 * @param minProb Minimum probability threshold
	 * @param minQual Minimum quality threshold
	 * @param allowZeroSizeSketch Whether to allow empty sketches
	 * @return Generated sketch from all processing threads
	 */
	public Sketch processReadsMT(FileFormat ffin1, FileFormat ffin2, int maxThreads, long reads, int mode, float samplerate,
			float minEntropy, float minProb, byte minQual, boolean allowZeroSizeSketch){
		assert(mode==ONE_SKETCH || mode==PER_FILE);
		final boolean compressed=ffin1.compressed();
		
		maxThreads=Tools.mid(1, maxThreads, Shared.threads());
		
		//Create a read input stream
		final Streamer cris;
		String simpleName;
		{
			simpleName=ffin1.simpleName();
			cris=StreamerFactory.getReadInputStream(reads, true, ffin1, ffin2, null, null, 1);
			if(samplerate!=1){cris.setSampleRate(samplerate, sampleseed);}
			cris.start(); //Start the stream
//			if(verbose){outstream.println("Started cris");}
		}
		
		final int threads=(int)Tools.min(maxThreads,
				2*(4)*(ffin2==null ? 4 : 8)*(mergePairs ? 3 : minEntropy>0 ? 2 : 1));
		
		if(verbose2 || true){System.err.println("Starting "+threads+" load threads.");}
		ArrayList<LoadThread2> list=new ArrayList<LoadThread2>(threads);
		for(int i=0; i<threads; i++){
			list.add(new LoadThread2(cris, minEntropy, minProb, minQual));
			list.get(i).start();
		}
		
		ArrayList<SketchHeap> heaps=new ArrayList<SketchHeap>(threads);
		
		for(LoadThread2 pt : list){

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
			
			if(pt.smm.heap!=null && pt.smm.heap.size()>0){
				heaps.add(pt.smm.heap);
			}
		}
		list.clear();
		ReadWrite.closeStream(cris);
		
		if(verbose2){System.err.println("Generating a sketch by combining thread output.");}
		Sketch sketch=toSketch(heaps, allowZeroSizeSketch);
		if(verbose2){System.err.println("Resulting sketch: "+((sketch==null) ? "null" : "length="+sketch.length()));}
		if(sketch!=null){sketch.setFname(simpleName);}
		return sketch;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Writing            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Writes sketches to multiple output files in round-robin fashion.
	 * @param sketches List of sketches to write
	 * @param ff Array of output file formats
	 * @return true if any errors occurred during writing
	 */
	public static boolean write(ArrayList<Sketch> sketches, FileFormat ff[]){
		final int len=ff.length;
		ByteStreamWriter tsw[]=new ByteStreamWriter[len];
		for(int i=0; i<len; i++){
			tsw[i]=new ByteStreamWriter(ff[i]);
			tsw[i].start();
		}
		boolean error=false;
		for(int i=0; i<sketches.size(); i++){
			write(sketches.get(i), tsw[i%len], new ByteBuilder());
		}
		for(int i=0; i<len; i++){
			error|=tsw[i].poisonAndWait();
		}
		return error;
	}
	
	/**
	 * Writes all sketches to a single output file.
	 * @param sketches List of sketches to write
	 * @param ff Output file format
	 * @return true if errors occurred during writing
	 */
	public static boolean write(ArrayList<Sketch> sketches, FileFormat ff){
		final ByteStreamWriter tsw=new ByteStreamWriter(ff);
		final ByteBuilder bb=new ByteBuilder();
		tsw.start();
		for(Sketch sketch : sketches){
			write(sketch, tsw, bb);
		}
		return tsw.poisonAndWait();
	}
	
	/**
	 * Writes a single sketch to an output file.
	 * @param sketch Sketch to write
	 * @param ff Output file format
	 * @return true if errors occurred during writing
	 */
	public static boolean write(Sketch sketch, FileFormat ff){
//		System.err.println(ff.name()+", "+new File(ff.name()).exists());
		ByteStreamWriter tsw=new ByteStreamWriter(ff);
//		assert(false) : new File(ff.name()).exists();
		tsw.start();
		write(sketch, tsw, null);
		return tsw.poisonAndWait();
	}
	
	/**
	 * Core sketch writing method using byte stream writer.
	 * @param sketch Sketch to write
	 * @param tsw Target byte stream writer
	 * @param bb Byte builder for temporary storage (may be null)
	 */
	public static void write(Sketch sketch, ByteStreamWriter tsw, ByteBuilder bb){
		if(bb==null){bb=new ByteBuilder();}
		else{bb.clear();}
		sketch.toBytes(bb);
		tsw.print(bb);
	}
		
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
//	final EntropyTracker eTracker;
	/** Target sketch size (maximum number of k-mers to retain) */
	final int stTargetSketchSize;
	/** Minimum k-mer occurrence count threshold for inclusion in sketch */
	public final int minKeyOccuranceCount;
	/** Force kmer counts to be tracked. */
	public final boolean trackCounts;
	/** Merge reads before processing kmers. */
	public final boolean mergePairs;
	/** Whether to include reverse complement k-mers in sketches */
	public final boolean rcomp;
	
	/** Buffer length for buffered sketch file reading */
	public static int BUFLEN=16384;
	/** Whether to use buffered reader for sketch file loading */
	public static boolean BUFFERED_READER=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
//	public static boolean verbose=false;
	
}
