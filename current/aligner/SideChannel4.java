package aligner;

import java.util.ArrayList;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBDuk;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.SamLine;
import stream.Writer;
import stream.WriterFactory;

/**
 * Side channel mapper for performing lightweight alignment of reads to reference sequences.
 * Uses one or two k-mer-based indices for fast mapping.
 * * Updated to use stream.Writer instead of ConcurrentReadOutputStream to match BBDukS architecture.
 * * @author Brian Bushnell
 */
public class SideChannel4 {
	
	/*--------------------------------------------------------------*/
	/*----------------          Constructor         ----------------*/
	/*--------------------------------------------------------------*/
	
	public SideChannel4(String ref_, String out_, String outu_, int k1_, float minid1, 
			int midMaskLen1, boolean overwrite_, boolean ordered_) {
		this(ref_, out_, outu_, k1_, -1, minid1, 1, midMaskLen1, 0, overwrite_, ordered_);
	}
	
	public SideChannel4(String ref_, String out_, String outu_, int k1_, int k2_, float minid1, float minid2, 
			int midMaskLen1, int midMaskLen2, boolean overwrite_, boolean ordered_) {
		Timer t=new Timer();
		ref=fixRefPath(ref_);
		out=out_;
		outu=outu_;
		k1=Tools.max(k1_, k2_);
		k2=Tools.min(k1_, k2_);
		minIdentity1=fixID(minid1);
		minIdentity2=fixID(minid2);
		overwrite=overwrite_;
		ordered=false;//ordered_;
		assert(k1>0);

		ffout=FileFormat.testOutput(out, FileFormat.SAM, null, true, overwrite, false, ordered);
		ffoutu=FileFormat.testOutput(outu, FileFormat.FASTQ, null, true, overwrite, false, ordered);
		samOut=((ffout!=null && ffout.samOrBam()) || ((ffoutu!=null && ffoutu.samOrBam())));
		
		final Read r=MicroIndex3.loadRef(ref, samOut);
		// Capture refName to avoid Global lookup crashes if Data.* is not initialized
		refName=(r!=null ? r.id : null);
		refNameB=(refName!=null ? refName.getBytes() : null);
		
		index1=new MicroIndex3(k1, midMaskLen1, r);
		index2=(k2<1 ? null : new MicroIndex3(k2, midMaskLen2, r));
		mapper1=new MicroAligner3(index1, minIdentity1, true);
		mapper2=(k2<1 ? null : new MicroAligner3(index2, minIdentity2, true));

		if(samOut) {ReadStreamByteWriter.USE_ATTACHED_SAMLINE=true;}
		
		final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
		
		// Use WriterFactory instead of ConcurrentReadOutputStream
		if(ffout!=null) {
			ros=WriterFactory.getStream(ffout, null, null, null, buff, null, ordered, -1);
//			ros.start();
		}else {ros=null;}
		
		if(ffoutu!=null) {
			rosu=WriterFactory.getStream(ffoutu, null, null, null, buff, null, ordered, -1);
//			rosu.start();
		}else {rosu=null;}
		
		t.stop("Created side channel"+(out==null ? "" : (" for "+out))+": ");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean map(Read r1, Read r2) {
		return map(r1, r2, mapper1, mapper2);
	}
	
	public boolean map(Read r1, Read r2, MicroAligner3 mapper1, MicroAligner3 mapper2) {
		float id1=mapper1.map(r1);
		float id2=mapper1.map(r2);
		if(id1+id2<=0) {return false;}//Common case
		
		if(r2!=null) {
			if(mapper2!=null) {
				if(r1.mapped() && !r2.mapped()) {id2=mapper2.map(r2);}
				else if(r2.mapped() && !r1.mapped()) {id1=mapper2.map(r1);}
			}
			boolean properPair=(r1.mapped() && r2.mapped() && r1.chrom==r2.chrom && 
					r1.strand()!=r2.strand() && Tools.absdif(r1.start, r2.start)<=1000);
			r1.setPaired(properPair);
			r2.setPaired(properPair);
		}
		
		if(!r1.mapped()) {id1=0;}
		if(r2==null || !r2.mapped()) {id2=0;}
		long idsum=(long)((id1+id2)*10000);
		if(idsum<=0) {return false;}
		
		return true;
	}
	
	public String stats(long readsIn, long basesIn) {
		long ro, bo, idsum, rm;
		ro=readsOut; bo=basesOut; idsum=identitySum; rm=readsMapped;
		String s=("Aligned reads:          \t"+ro+" reads ("+BBDuk.toPercent(ro, readsIn)+") \t"+
				+bo+" bases ("+BBDuk.toPercent(bo, basesIn)+") \tavgID="+Tools.format("%.4f", idsum/(100.0*rm)));
		return s;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             I/O              ----------------*/
	/*--------------------------------------------------------------*/
	
	public void start() {
		if(ros!=null) {ros.start();}
		if(rosu!=null) {rosu.start();}
	}
	
	public boolean shutdown() {
		// Updated to close Writer streams
		errorState=ReadWrite.closeOutputStreams(ros, rosu)|errorState;
		return errorState;
	}

	public int writeToMapped(ArrayList<Read> reads, long num) {
		// IMPORTANT: Removed the optimization that skips empty lists. 
		// Writer needs to see the ID increment even for empty batches.
		if(reads==null || (ros==null && !TRACK_STATS)) {
			return 0;
		}
		
		int rm=0, ro=0, bo=0;
		double idsum=0;
		for(Read r1 : reads) {
			Read r2=r1.mate;
			ro+=r1.pairCount();
			bo+=r1.pairLength();
			rm+=r1.pairMappedCount();

			if(TRACK_STATS) {
				if(r1.mapped()) {idsum+=r1.identity();}
				if(r2!=null && r2.mapped()) {idsum+=r2.identity();}
			}
		}

		if(TRACK_STATS && ro>0) {
			synchronized(this) {
				readsMapped+=rm;
				readsOut+=ro;
				basesOut+=bo;
				identitySum+=(long)(idsum*10000);
			}
		}
		if(ros==null) {return ro;}
		
		final boolean makeSamLine=(samOut && ReadStreamByteWriter.USE_ATTACHED_SAMLINE);
		if(makeSamLine) {
			for(Read r1 : reads) {
				Read r2=r1.mate;
				r1.samline=(r1.samline!=null ? r1.samline : new SamLine(r1, 0));
				
				// Added safety: Set RNAME from local refName to avoid Global lookup
				if(r1.mapped() && refName!=null) {
					r1.samline.setRname(refNameB);
				}
				
				if(r2!=null) {
					r2.samline=(r2.samline!=null ? r2.samline : new SamLine(r2, 1));
					if(r2.mapped() && refName!=null) {
						r2.samline.setRname(refNameB);
					}
				}
			}
		}
		// Using Writer.add instead of ConcurrentReadOutputStream.add
		ros.add(reads, num);
		return ro;
	}

	public int writeByStatus(ArrayList<Read> reads, long num) {
		int ro=writeMappedOnly(reads, num);
		writeUnmappedOnly(reads, num);
		return ro;
	}

	private int writeMappedOnly(ArrayList<Read> reads, long num) {
		// Always pass through to Writer to maintain ID sequence
		if(reads==null) {return 0;}
		
		if(ros==null) {
			// If no output stream, we still need to track stats, but we can't write.
			// However, if ros IS null, we don't need to worry about the ID sequence for it.
			// But we DO need to check if TRACK_STATS is on.
			if(!TRACK_STATS) return 0;
		}
		
		int listSize=0;
		int rm=0, ro=0, bo=0, rou=0, bou=0;
		double idsum=0;
		for(Read r1 : reads) {
			boolean mapped=r1.eitherMapped();
			if(mapped) {
				Read r2=r1.mate;
				listSize++;
				ro+=r1.pairCount();
				bo+=r1.pairLength();
				rm+=r1.pairMappedCount();

				if(TRACK_STATS) {
					if(r1.mapped()) {idsum+=r1.identity();}
					if(r2!=null && r2.mapped()) {idsum+=r2.identity();}
				}
			}else{
				rou+=r1.pairCount();
				bou+=r1.pairLength();
			}
		}

		if(TRACK_STATS && ro>0) {
			synchronized(this) {
				readsMapped+=rm;
				readsOut+=ro;
				basesOut+=bo;
				identitySum+=(long)(idsum*10000);
			}
		}

		if(ros==null) {return ro;}

		ArrayList<Read> list=new ArrayList<Read>(Tools.max(1, listSize));
		final boolean makeSamline=(samOut && ReadStreamByteWriter.USE_ATTACHED_SAMLINE);
		for(Read r1 : reads) {
			boolean mapped=r1.eitherMapped();
			if(mapped && list!=null) {
				Read r2=r1.mate;
				if(makeSamline) {
					r1.samline=(r1.samline!=null ? r1.samline : new SamLine(r1, 0));
					if(r1.mapped() && refName!=null) r1.samline.setRname(refNameB);
					
					if(r2!=null) {
						r2.samline=(r2.samline!=null ? r2.samline : new SamLine(r2, 1));
						if(r2.mapped() && refName!=null) r2.samline.setRname(refNameB);
					}
				}
				list.add(r1);
			}
		}
		ros.add(list, num);
		return ro;
	}

	private int writeUnmappedOnly(ArrayList<Read> reads, long num) {
		if(reads==null || rosu==null) {
			return 0;
		}

		ArrayList<Read> list=new ArrayList<Read>(8);
		int rou=0;
		for(Read r1 : reads) {
			if(!r1.eitherMapped()) {
				rou+=r1.pairCount();
				list.add(r1);
			}
		}
		rosu.add(list, num);
		return rou;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int[] parseK(String arg, String a, String b) {
		int[] ret=new int[2];
		String[] terms=b.split(",");
		for(int i=0; i<terms.length; i++) {
			ret[i]=Integer.parseInt(terms[i]);
		}
		return ret;
	}
	
	static float fixID(float id) {
		if(id>1) {id=id/100;}
		assert(id<=1);
		return id;
	}
	
	static String fixRefPath(String refPath) {
		if(refPath==null || Tools.isReadableFile(refPath)) {return refPath;}
		if("phix".equalsIgnoreCase(refPath)){return Data.findPath("?phix2.fa.gz");}
		return refPath;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public final MicroIndex3 index1;
	public final MicroIndex3 index2;
	public final MicroAligner3 mapper1;
	public final MicroAligner3 mapper2;
	
	public final int k1;
	public final int k2;
	public final float minIdentity1;
	public final float minIdentity2;
	
	public final String ref;
	public final String out;
	public final String outu;
	public final boolean samOut;
	public final FileFormat ffout;
	public final FileFormat ffoutu;
	
	// Changed from ConcurrentReadOutputStream to Writer
	private final Writer ros;
	private final Writer rosu;

	private final String refName;
	private final byte[] refNameB;
	
	public long readsMapped=0;
	public long readsOut=0;
	public long basesOut=0;
	public long identitySum=0;//x100%; 0-10000 scale. 
	
	public final boolean overwrite;
	public final boolean ordered;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean TRACK_STATS=true;

}