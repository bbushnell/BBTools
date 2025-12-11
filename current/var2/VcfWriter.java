package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import gff.GffLine;
import ml.CellNet;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Multithreaded writer for variant data in VCF, VAR, or GFF formats.
 * Uses producer-consumer pattern with multiple threads formatting variants while a single ByteStreamWriter handles ordered output to maintain proper file structure.
 * Supports comprehensive filtering, statistical metadata inclusion, and format-specific optimizations for each output type.
 * @author Brian Bushnell
 * @contributor Isla
 */
public class VcfWriter {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	public VcfWriter(VarMap varMap_, VarFilter filter_, long reads_,
			long pairs_, long properPairs_, long bases_, String ref_,
			boolean trimWhitespace_, String sampleName_){
		
		varMap=varMap_;
		filter=filter_;
		trimWhitespace=trimWhitespace_;
		
		reads=reads_;
		pairs=pairs_;
		properPairs=properPairs_;
		bases=bases_;
		ref=ref_;
		sampleName=sampleName_;
		
		ploidy=varMap.ploidy;
		properPairRate=varMap.properPairRate;
		pairedInSequencingRate=varMap.pairedInSequencingRate;
		totalQualityAvg=varMap.totalQualityAvg;
		totalMapqAvg=varMap.totalMapqAvg;
		readLengthAvg=varMap.readLengthAvg;
		scafMap=varMap.scafMap;
		
		threads=Tools.max(1, Shared.threads());
		
		array=varMap.toArray(true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean writeVcf(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeVcfFile(fname);
		return vw.errorState;
	}
	
	public static boolean writeVar(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeVarFile(fname);
		return vw.errorState;
	}
	
	public static boolean writeGff(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeGffFile(fname);
		return vw.errorState;
	}

	/** Writes variants to VCF format file.
	 * @param fname Output filename */
	public void writeVcfFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.VCF, "vcf", true, true, false, true);
		writeFile(ff, VCFMODE);
	}

	/** Writes variants to VAR format file.
	 * @param fname Output filename */
	public void writeVarFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.VAR, "var", true, true, false, true);
		writeFile(ff, VARMODE);
	}

	/** Writes variants to GFF format file.
	 * @param fname Output filename */
	public void writeGffFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.GFF, "gff", true, true, false, true);
		writeFile(ff, GFFMODE);
	}

	/** Writes variants to VCF format using FileFormat.
	 * @param ff Output file format */
	public void writeVcfFile(final FileFormat ff){
		assert(ff.vcf());
		writeFile(ff, VCFMODE);
	}

	/** Writes variants to VAR format using FileFormat.
	 * @param ff Output file format */
	public void writeVarFile(final FileFormat ff){
		assert(ff.var()) : "Incorrect file extension: "+ff;
		writeFile(ff, VARMODE);
	}

	/** Writes variants to GFF format using FileFormat.
	 * @param ff Output file format */
	public void writeGffFile(final FileFormat ff){
		assert(ff.gff());
		writeFile(ff, GFFMODE);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Core file writing method using multithreaded processing.
	 * Coordinates header writing, thread management, and output ordering.
	 * @param ff Output file format
	 * @param mode Output format mode (VCF, VAR, or GFF)
	 */
	private void writeFile(final FileFormat ff, final int mode){
		assert(ff.ordered());
		final ArrayBlockingQueue<ListNum<Var>> inq=new ArrayBlockingQueue<ListNum<Var>>(threads+1);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		
		writeHeader(bsw, mode);
		
		ArrayList<ProcessThread> alpt=spawnThreads(bsw, inq, mode);
		
		makeLists(inq);
		
		waitForFinish(alpt);
		
		errorState=bsw.poisonAndWait()|errorState;
	}
	
	private void writeHeader(ByteStreamWriter bsw, int mode){
		ByteBuilder bb=new ByteBuilder(1000);
		if(mode==VCFMODE){
			bb.append(VarHelper.toVcfHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					ploidy, reads, pairs, properPairs, bases, ref, scafMap, sampleName, trimWhitespace)).append('\n');
		}else if(mode==VARMODE){
			bb.append(VarHelper.toVarHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					ploidy, reads, pairs, properPairs, bases, ref)).append('\n');
		}else if(mode==GFFMODE){
			bb.append(GffLine.toHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					ploidy, reads, pairs, properPairs, bases, ref)).append('\n');
		}else{assert(false);}
		bsw.add(bb, 0);
	}

	/**
	 * Spawns worker threads for variant formatting.
	 * Each thread formats variants from the queue into the appropriate output format.
	 *
	 * @param bsw ByteStreamWriter for ordered output
	 * @param inq Input queue for variant batches
	 * @param mode Output format mode
	 * @return List of spawned ProcessThread objects
	 */
	private ArrayList<ProcessThread> spawnThreads(ByteStreamWriter bsw, final ArrayBlockingQueue<ListNum<Var>> inq, int mode){
		
		//Do anything necessary prior to processing
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, bsw, inq, mode));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		return alpt;
	}
	
	private void waitForFinish(ArrayList<ProcessThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(ProcessThread pt : alpt){
			if(verbose){outstream.println("Waiting for thread "+pt.tid);}
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}

			allSuccess&=pt.success;
		}
		
		//Track whether any threads failed
		if(!allSuccess){errorState=true;}
	}
	
	/**
	 * Creates batches of variants and queues them for processing.
	 * Distributes variants among worker threads while maintaining order.
	 * @param inq Queue for variant batches
	 */
	void makeLists(ArrayBlockingQueue<ListNum<Var>> inq){
		ArrayList<Var> list=new ArrayList<Var>(LIST_SIZE);
		long nextJobID=1;
		for(Var v : array){
			list.add(v);
			if(list.size()>=LIST_SIZE){
				putVars(new ListNum<Var>(list, nextJobID), inq);
				nextJobID++;
				list=new ArrayList<Var>(LIST_SIZE);
			}
		}
		if(list.size()>0){
			putVars(new ListNum<Var>(list, nextJobID), inq);
			nextJobID++;
			list=null;
		}
		if(verbose){outstream.println("tid "+0+" done making var lists.");}
		putVars(POISON_VARS, inq);
		if(verbose){outstream.println("tid "+0+" done poisoning.");}
	}
	
	final void putVars(ListNum<Var> list, final ArrayBlockingQueue<ListNum<Var>> inq){
		if(verbose){outstream.println("tid "+0+" putting vlist "+list.id+", size "+list.size());}
		while(list!=null){
			try {
				inq.put(list);
				list=null;
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		if(verbose){outstream.println("tid "+0+" done putting vlist");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread that formats variants into output format.
	 * Takes variant batches from queue, applies filtering, and formats
	 * each variant according to the specified output mode.
	 */
	private class ProcessThread extends Thread{
		
		public ProcessThread(int tid_, ByteStreamWriter bsw_, ArrayBlockingQueue<ListNum<Var>> inq_, int mode_){
			tid=tid_;
			bsw=bsw_;
			inq=inq_;
			mode=mode_;
		}
		
		/** Main thread execution - processes variant batches. */
		@Override
		/** Main thread execution - processes variant batches */
		public void run(){
			makeBytes();
		}
		
		/**
		 * Core variant processing loop.
		 * Takes variant batches, applies filtering, formats according to mode,
		 * and sends formatted output to ByteStreamWriter.
		 */
		void makeBytes(){
			if(verbose){outstream.println("tid "+tid+" started makeBytes.");}
			
			ListNum<Var> list=takeVars();
			final ByteBuilder bb=new ByteBuilder();
			while(list!=null && list.size()>0 && list!=POISON_VARS){
				
				for(Var v : list){
					if(v.forced() || filter==null || !filter.failNearby || v.nearbyVarCount<=filter.maxNearbyCount){
						if(mode==VCFMODE){
							v.toVCF(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap, filter, net, trimWhitespace);
						}else if(mode==VARMODE){
							v.toText(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap, net);
						}else if(mode==GFFMODE){
							GffLine.toText(bb, v, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap, net);
						}
						bb.nl();
					}
				}
				
				bsw.add(new ByteBuilder(bb.toBytes()), list.id);
				list=takeVars();
				
				bb.clear();
				bb.shrinkTo(SHRINK_SIZE);
			}
			if(list==POISON_VARS){
				putVars(POISON_VARS, inq);
			}
			if(verbose){outstream.println("tid "+tid+" done making bytes.");}
		}
		
		final ListNum<Var> takeVars(){
			if(verbose){outstream.println("tid "+tid+" taking vlist");}
			ListNum<Var> list=null;
			while(list==null){
				try {
					list=inq.take();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" took vlist "+list.id+", size "+list.size());}
			return list;
		}
		
		final ArrayBlockingQueue<ListNum<Var>> inq;
		final int tid;
		final ByteStreamWriter bsw;
		CellNet net;
		final int mode;
		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	final Var[] array;
	final int ploidy;
	final double properPairRate;
	final double pairedInSequencingRate;
	final double totalQualityAvg;
	final double totalMapqAvg;
	final double readLengthAvg;
	final ScafMap scafMap;
	final VarMap varMap;
	final VarFilter filter;
	final boolean trimWhitespace;
	
	final String sampleName;
	final long reads;
	final long pairs;
	final long properPairs;
	final long bases;
	final String ref;
	
	final int threads;
	
	boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int VARMODE=0, VCFMODE=1, GFFMODE=2;
	private static final ListNum<Var> POISON_VARS=new ListNum<Var>(null, Long.MAX_VALUE, true, false);
	
	private static boolean verbose=false;

	private static final int LIST_SIZE=200;
	private static final int SHRINK_SIZE=1000*LIST_SIZE;
	
	protected static PrintStream outstream=System.err;
}