package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import jgi.TranslateSixFrames;
import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import prok.ProkObject;
import shared.Parse;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import tax.TaxNode;
import tax.TaxTree;
import tracker.EntropyTracker;

/**
 * Creates MinHashSketches rapidly.
 * 
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public class SketchMakerMini extends SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public SketchMakerMini(SketchTool tool_, int mode_, DisplayParams params){
		this(tool_, mode_, params.minEntropy, params.minProb, params.minQual);
	}
	
	/**
	 * Constructor.
	 */
	public SketchMakerMini(SketchTool tool_, int mode_, float minEntropy_, float minProb_, byte minQual_){
		
		tool=tool_;
		mode=mode_;
		minProb=minProb_;
		minQual=minQual_;
		rcomp=tool.rcomp;
		
		aminoShift=AminoAcid.AMINO_SHIFT;
		if(!aminoOrTranslate()){
			shift=2*k;
			shift2=shift-2;
			mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		}else{
			shift=aminoShift*k;
			shift2=shift-aminoShift;
			mask=(shift>63 ? -1L : ~((-1L)<<shift));
		}
		
		if(AUTOSIZE && (mode==ONE_SKETCH || mode==PER_FILE)){
			heap=new SketchHeap(Tools.max(tool.stTargetSketchSize, (int)(80000*Tools.mid(1, AUTOSIZE_FACTOR, 32))), tool.minKeyOccuranceCount, tool.trackCounts);
		}else if(AUTOSIZE_LINEAR && (mode==ONE_SKETCH || mode==PER_FILE)){
			heap=new SketchHeap(Tools.max(tool.stTargetSketchSize, (int)(10000000*Tools.mid(0.1, 2*AUTOSIZE_LINEAR_DENSITY, 0.00001))), 
					tool.minKeyOccuranceCount, tool.trackCounts);
		}else{
			heap=new SketchHeap(tool.stTargetSketchSize, tool.minKeyOccuranceCount, tool.trackCounts);
		}
		
		if(minEntropy_>0){
			eTracker=new EntropyTracker(entropyK, k, (amino || translate), minEntropy_, true);
		}else{
			eTracker=null;
		}
		
		if(translate || processSSU){
			gCaller=CallGenes.makeGeneCaller(pgm);
		}else{
			gCaller=null;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public ArrayList<Sketch> toSketches(final String fname, float samplerate, long reads){
		final FileFormat ffin1, ffin2;
		if(fname.indexOf('#')>=0 && FileFormat.isFastqExt(ReadWrite.rawExtension(fname)) && !new File(fname).exists()){
			ffin1=FileFormat.testInput(fname.replaceFirst("#", "1"), FileFormat.FASTQ, null, true, true);
			ffin2=FileFormat.testInput(fname.replaceFirst("#", "2"), FileFormat.FASTQ, null, true, true);
		}else{
			ffin1=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
			ffin2=null;
		}
		return toSketches(ffin1, ffin2, samplerate, reads);
	}

	/** Create read streams and process all data */
	public ArrayList<Sketch> toSketches(FileFormat ffin1, FileFormat ffin2, float samplerate, long reads){
		heap.clear(false); //123
		final String simpleName;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			simpleName=ffin1.simpleName();
			heap.setFname(simpleName);
			cris=ConcurrentReadInputStream.getReadInputStream(reads, true, ffin1, ffin2, null, null);
			if(samplerate!=1){cris.setSampleRate(samplerate, sampleseed);}
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		if(mode==ONE_SKETCH || mode==PER_FILE){
			 if(heap.name0()==null){heap.setName0(simpleName);}
		}
		ArrayList<Sketch> sketches=processInner(cris);
		
		errorState|=ReadWrite.closeStream(cris);
		sketchesMade++;
		return sketches;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Iterate through the reads */
	ArrayList<Sketch> processInner(ConcurrentReadInputStream cris){
		assert(heap.size()==0);
		ArrayList<Sketch> sketches=new ArrayList<Sketch>(mode==ONE_SKETCH || mode==PER_FILE ? 1 : 8);
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();
		//Grab the actual read list from the ListNum
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//As long as there is a nonempty read list...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;
				
				processReadPair(r1, r2);
				if(mode!=ONE_SKETCH && mode!=PER_FILE){
					if(heap!=null && heap.size()>0 && heap.maxLen()>=Tools.max(1, minSketchSize)){
						int size=heap.size();
						Sketch sketch=new Sketch(heap, false, tool.trackCounts, null);
						assert(sketch.keys.length>0) : sketch.keys.length+", "+size;
						sketch.loadSSU();
						sketches.add(sketch);
						sketchesMade++;
					}
					if(heap!=null){heap.clear(false);}
				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
			//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

			//Fetch a new list
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		
		if(mode==ONE_SKETCH || mode==PER_FILE){
			Sketch sketch=new Sketch(heap, false, tool.trackCounts, null);
			sketch.loadSSU();
			sketches.add(sketch);
			sketchesMade++;
		}
		heap.clear(true);
		return sketches;
	}

	void processReadPair(Read r1, Read r2){
		//Track the initial length for statistics
		final int initialLength1=r1.length();
		final int initialLength2=r1.mateLength();

		//Increment counters
		readsProcessed+=r1.pairCount();
		basesProcessed+=initialLength1+initialLength2;
		
		if(mode!=ONE_SKETCH && mode!=PER_FILE){
			int expectedSize=toSketchSize(initialLength1+initialLength2, -1, -1, targetSketchSize);
			if(heap==null || heap.capacity()<expectedSize){heap=new SketchHeap(expectedSize, tool.minKeyOccuranceCount, tool.trackCounts);}
		}
		
		if(tool.mergePairs && r2!=null){
			final int insert=BBMerge.findOverlapStrict(r1, r2, false);
			if(insert>0){
				heap.genomeSequences++;
				r2.reverseComplement();
				r1=r1.joinRead(insert);
				r2=null;
			}
		}
		
		processRead(r1);
		if(r2!=null){processRead(r2);}

		if(heap.name0()==null){
			heap.setName0(r1.id);
		}
		
		TaxNode tn=null;
		if(heap.taxID<0 && r1.length()>800){
			if(taxtree!=null){ 
				try {
					tn=taxtree.parseNodeFromHeader(r1.id, true);
				} catch (Throwable e) {}
				if(tn!=null){
					heap.taxID=tn.id;
					if(heap.taxName()==null){
						heap.setTaxName(tn.name);
					}
				}
//				System.err.println("A) "+heap.taxID+r1.id);
			}else{
				heap.taxID=TaxTree.parseHeaderStatic(r1.id);
//				System.err.println("B) "+heap.taxID+r1.id);
			}
		}
		assert(heap.taxID<0 || heap.taxName()!=null || taxtree==null) : heap.taxID+", "+heap.taxName()+", "+heap.name()+", "+tn;
	}

	public void processRead(final Read r){
		if(amino){
			processReadAmino(r);
		}else if(translate){
			processReadTranslated(r);
		}else{
			processReadNucleotide(r);
		}
	}
	
	public void processReadTranslated(final Read r){
		assert(!r.aminoacid());
		final ArrayList<Read> prots;
		if(sixframes){
			if(processSSU && heap.r16S()==null && r.length()>=min_SSU_len && !useSSUMapOnly && !heap.isEukaryote()){
				Orf orf=gCaller.makeRna(r.id, r.bases, ProkObject.r16S);//TODO: allow 18S also
				if(orf!=null && orf.length()>=min_SSU_len){
					assert(orf.is16S());
					if(orf.is16S() && orf.length()>=heap.r16SLen()){heap.set16S(CallGenes.fetch(orf, r).bases);}
				}
				//TODO: Add 18S.
			}
			prots=TranslateSixFrames.toFrames(r, true, false, 6);
		}else{
			ArrayList<Orf> list;
			list=gCaller.callGenes(r);
			prots=CallGenes.translate(r, list);
			if(processSSU && heap.r16S()==null && r.length()>=min_SSU_len && !useSSUMapOnly && !heap.isEukaryote()){
				for(Orf orf : list){
					if(orf.is16S() && orf.length()>=min_SSU_len && orf.length()>=heap.r16SLen()){
						heap.set16S(CallGenes.fetch(orf, r).bases);
						break;
					}
				}
			}
		}
		if(prots!=null){
			for(Read p : prots){
				processReadAmino(p);
			}
		}
	}
	
	public void processReadNucleotide(final Read r){
		if(processSSU && heap.r16S()==null && r.length()>=min_SSU_len && !useSSUMapOnly && !heap.isEukaryote()){
			Orf orf=gCaller.makeRna(r.id, r.bases, ProkObject.r16S);//TODO: 18S
			if(orf!=null && orf.length()>=min_SSU_len){
				assert(orf.start>=0 && orf.stop<r.length()) : r.length()+"\n"+orf;
				assert(orf.is16S());
				if(orf.is16S() && orf.length()>=heap.r16SLen()){heap.set16S(CallGenes.fetch(orf, r).bases);}
			}
			//TODO: Add 18S.
		}
		
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final long[] baseCounts=heap.baseCounts(true);
		long kmer=0;
		long rkmer=0;
		int len=0;
		assert(!r.aminoacid());
		
		final boolean noBlacklist=!(Blacklist.exists() || Whitelist.exists());
		final long min=minHashValue;
		heap.genomeSizeBases+=r.length();
		heap.genomeSequences++;
		if(eTracker!=null){eTracker.clear();}
		
//		assert(false) : minProb+", "+minQual+", "+(quals==null);
		
		if(quals==null || (minProb<=0 && minQual<2)){
//			System.err.println("A");
			for(int i=0; i<bases.length; i++){
//				System.err.println("B: len="+len);
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(eTracker!=null){eTracker.add(b);}
				if(x<0){
					len=0;
					rkmer=0;
				}else{
					len++;
					baseCounts[(int)x]++;
				}
				
//				System.err.println("\n"+AminoAcid.kmerToString(kmer, k)+"\n"+AminoAcid.kmerToString(rkmer, k)+"\n"
//				+AminoAcid.kmerToString(AminoAcid.reverseComplementBinaryFast(rkmer, k), k)+"\n"
//				+len+", "+(char)b+", "+x+", "+x2+"\n");
				
				if(len>=k){
					kmersProcessed++;
					heap.genomeSizeKmers++;
//					heap.probSum++; //Note really necessary for fasta data
					if(eTracker==null || eTracker.passes()){
//						System.err.println("Pass.\t"+eTracker.calcEntropy()+"\t"+eTracker.basesToString());
						
//						assert(kmer==AminoAcid.reverseComplementBinaryFast(rkmer, k)) : 
//							"\n"+AminoAcid.kmerToString(kmer, k)+"\n"+AminoAcid.kmerToString(rkmer, k)+"\n"
//							+AminoAcid.kmerToString(AminoAcid.reverseComplementBinaryFast(rkmer, k), k)+"\n"
//							+len+", "+(char)b+", "+x+", "+x2;
						
						final long hashcode=hash(kmer, rkmer, rcomp);
						//				System.err.println(kmer+"\t"+rkmer+"\t"+z+"\t"+hash);
						if(hashcode>min){
							if(noBlacklist){
								heap.add(hashcode);
							}else{
								heap.checkAndAdd(hashcode);
							}
						}
					}else{
//						System.err.println("Fail.\t"+eTracker.calcEntropy()+"\t"+eTracker.basesToString()+"\n"+r.toFastq()+"\n"+eTracker);
//						assert(false);
					}
				}
			}
		}else{
			int zeroQualityKmers=0;
			int positiveQualityKmers=0;
			
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];

				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(eTracker!=null){eTracker.add(b);}

				final byte q=quals[i];
				{//Quality-related stuff
					assert(q>=0) : Arrays.toString(quals)+"\n"+minProb+", "+minQual;
//					if(x>=0){
//						if(q>0){
//							positiveQualityBases++;
//						}else{
//							zeroQualityBases++;
//						}
//					}
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(len>k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
					if(x<0 || q<minQual){
						len=0;
						kmer=rkmer=0;
						prob=1;
					}else{
						len++;
						baseCounts[(int)x]++;
					}
				}
				
				if(len>=k){
					kmersProcessed++;
					if(prob>=minProb){
						heap.genomeSizeKmers++;
						heap.probSum+=prob;
						if(eTracker==null || eTracker.passes()){
							final long hashcode=hash(kmer, rkmer, rcomp);
//							System.err.println(kmer+"\t"+rkmer+"\t"+z+"\t"+hash);
							if(hashcode>min){
								if(noBlacklist){
									heap.add(hashcode);
								}else{
									heap.checkAndAdd(hashcode);
								}
							}
						}else{
//							System.err.println("Fail.\t"+eTracker.calcEntropy()+"\t"+eTracker.basesToString());
						}
						positiveQualityKmers++;
					}else if(q<=2){
						zeroQualityKmers++;
					}
				}
				
				//This version is slow but calculates depth better.
//				if(len>=k){
//					kmersProcessed++;
//					heap.genomeSizeKmers++;
//					final long hash=hash(kmer, rkmer, rcomp);
//					//				System.err.println(kmer+"\t"+rkmer+"\t"+z+"\t"+hash);
//					if(hash>min){
//						if(prob>=minProb || (!heap.setMode && heap.contains(hash))){
//							if(noBlacklist){
//								heap.add(hash);
//							}else{
//								heap.checkAndAdd(hash);
//							}
//						}
//					}
//				}
			}
			if(minProb>0 && zeroQualityKmers>100 && positiveQualityKmers==0){
				if(looksLikePacBio(r)){
					synchronized(this){
						minProb=0;
					}
					processReadNucleotide(r);
				}
			}
		}
//		assert(false);
	}
	
	boolean looksLikePacBio(Read r){
		if(r.length()<302 || r.mate!=null){return false;}
		if(r.quality==null){
			int x=Parse.parseZmw(r.id);
			return x>=0;
		}
		int positive=0;
		int zero=0;
		int ns=0;
		for(int i=0; i<r.bases.length; i++){
			byte b=r.bases[i];
			byte q=r.quality[i];
			if(b=='N'){ns++;}
			else if(q==0 || q==2){
				zero++;
			}else{
				positive++;
			}
		}
		return zero>=r.length()/2 && positive==0;
	}

	void processReadAmino(final Read r){
		final byte[] bases=r.bases;
		long kmer=0;
		int len=0;
		assert(r.aminoacid());
		
		final boolean noBlacklist=!(Blacklist.exists() || Whitelist.exists());
		final long min=minHashValue;
		heap.genomeSizeBases+=r.length();
		heap.genomeSequences++;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.acidToNumberNoStops[b];
			kmer=((kmer<<aminoShift)|x)&mask;
//			if(eTracker!=null){eTracker.add(b);}
			
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				kmersProcessed++;
				heap.genomeSizeKmers++;
//				if(eTracker==null || eTracker.passes()){
//					assert(false) : (eTracker==null)+", "+eTracker.cutoff()+", "+eTracker.calcEntropy()+", "+r;
					long hashcode=hash(kmer, kmer);
					if(hashcode>min){
						if(noBlacklist){
							heap.add(hashcode);
						}else{
							heap.checkAndAdd(hashcode);
						}
					}
//				}
			}
		}
	}

	void processReadAmino_old_no_entropy(final Read r){
		final byte[] bases=r.bases;
		long kmer=0;
		int len=0;
		assert(r.aminoacid());
		
		final boolean noBlacklist=!(Blacklist.exists() || Whitelist.exists());
		final long min=minHashValue;
		heap.genomeSizeBases+=r.length();
		heap.genomeSequences++;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.acidToNumberNoStops[b];
			kmer=((kmer<<aminoShift)|x)&mask;
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				kmersProcessed++;
				heap.genomeSizeKmers++;
				long hashcode=hash(kmer, kmer);
				if(hashcode>min){
					if(noBlacklist){
						heap.add(hashcode);
					}else{
						heap.checkAndAdd(hashcode);
					}
				}
			}
		}
	}
	
	public Sketch toSketch(int minCount){
		Sketch sketch=null;
		if(heap!=null && heap.size()>0){
			try {
				sketch=new Sketch(heap, false, tool.trackCounts, null, minCount);
			} catch (Throwable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			heap.clear(false);
		}
		return sketch;
	}
	
	public void add(SketchMakerMini smm){
		heap.add(smm.heap);
		readsProcessed+=smm.readsProcessed;
		basesProcessed+=smm.basesProcessed;
		kmersProcessed+=smm.kmersProcessed;
		sketchesMade+=smm.sketchesMade;
		pacBioDetected|=smm.pacBioDetected;
	}

	public boolean isEmpty() {
		return heap==null || heap.size()<1;
	}
	
	public float minEntropy() {
		// TODO Auto-generated method stub
		return eTracker==null ? -1 : eTracker.cutoff();
	}
	
	public SketchHeap heap() {return heap;}
	
//	public void clear() {
//		heap.clear();
//		readsProcessed=basesProcessed=kmersProcessed=sketchesMade=0;
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	

	/** True only if this thread has completed successfully */
	boolean success=false;

	SketchHeap heap;

	final int aminoShift;
	final int shift;
	final int shift2;
	final long mask;
	final EntropyTracker eTracker;
	final GeneCaller gCaller;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	public long readsProcessed=0;
	/** Number of bases processed */
	public long basesProcessed=0;
	/** Number of bases processed */
	public long kmersProcessed=0;
	/** Number of sketches started */
	public long sketchesMade=0;

	float minProb() {return minProb;}
	byte minQual() {return minQual;}
	public boolean pacBioDetected=false;
	private float minProb;
	private byte minQual;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final SketchTool tool;
	private final boolean rcomp;
	final int mode;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
