package aligner;

import java.util.ArrayList;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBDuk;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.SamLine;

public class SideChannel3 {
	
	/*--------------------------------------------------------------*/
	/*----------------          Constructor         ----------------*/
	/*--------------------------------------------------------------*/
	
	public SideChannel3(String ref_, String out_, String outu_, int k1_, float minid1, 
			int midMaskLen1, boolean overwrite_, boolean ordered_) {
		this(ref_, out_, outu_, k1_, -1, minid1, 1, midMaskLen1, 0, overwrite_, ordered_);
	}
	
	public SideChannel3(String ref_, String out_, String outu_, int k1_, int k2_, float minid1, float minid2, 
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
		ordered=ordered_;
		assert(k1>0);

		ffout=FileFormat.testOutput(out, FileFormat.SAM, null, true, overwrite, false, ordered);
		ffoutu=FileFormat.testOutput(outu, FileFormat.FASTQ, null, true, overwrite, false, ordered);
		samOut=((ffout!=null && ffout.samOrBam()) || ((ffoutu!=null && ffoutu.samOrBam())));
		final Read r=MicroIndex3.loadRef(ref, samOut);
		index1=new MicroIndex3(k1, midMaskLen1, r);
		index2=(k2<1 ? null : new MicroIndex3(k2, midMaskLen2, r));
		mapper1=new MicroAligner3(index1, minIdentity1, true);
		mapper2=(k2<1 ? null : new MicroAligner3(index2, minIdentity2, true));

		if(samOut) {ReadStreamByteWriter.USE_ATTACHED_SAMLINE=true;}
		final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
		if(ffout!=null) {
			cros=ConcurrentReadOutputStream.getStream(ffout, null, buff, null, false);
		}else {cros=null;}
		if(ffoutu!=null) {
			crosu=ConcurrentReadOutputStream.getStream(ffoutu, null, buff, null, false);
		}else {crosu=null;}
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
		if(cros!=null) {cros.start();}
		if(crosu!=null) {crosu.start();}
	}
	
	public boolean shutdown() {
		errorState=ReadWrite.closeOutputStreams(cros, crosu)|errorState;
		return errorState;
	}

	/** Everything in this list gets written */
	public int writeToMapped(ArrayList<Read> reads, long num) {
		if(reads==null || (reads.isEmpty() && !ordered) || (cros==null && !TRACK_STATS)) {
			assert(!ordered || reads!=null);
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
		if(cros==null) {return ro;}
		
		final boolean makeSamLine=(samOut && ReadStreamByteWriter.USE_ATTACHED_SAMLINE);
		if(makeSamLine) {
			for(Read r1 : reads) {
				Read r2=r1.mate;
				r1.samline=(r1.samline!=null ? r1.samline : new SamLine(r1, 0));
				if(r2!=null) {r2.samline=(r2.samline!=null ? r2.samline : new SamLine(r2, 1));}
			}
		}
		cros.add(reads, num);
		return ro;
	}

	/** Expects a list of mixed mapped and unmapped reads */
	public int writeByStatus(ArrayList<Read> reads, long num) {
		int ro=writeMappedOnly(reads, num);
		writeUnmappedOnly(reads, num);
		return ro;
	}

	/** Expects a list of mixed mapped and unmapped reads; 
	 * only writes mapped reads (plus mates) */
	private int writeMappedOnly(ArrayList<Read> reads, long num) {
		if(reads==null || (reads.isEmpty() && !ordered)) {
			assert(!ordered || reads!=null);
			return 0;
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

		if(cros==null || (!ordered && listSize<1)) {return ro;}

		ArrayList<Read> list=new ArrayList<Read>(Tools.max(1, listSize));
		final boolean makeSamline=(samOut && ReadStreamByteWriter.USE_ATTACHED_SAMLINE);
		for(Read r1 : reads) {
			boolean mapped=r1.eitherMapped();
			if(mapped && list!=null) {
				Read r2=r1.mate;
				if(makeSamline) {
					r1.samline=(r1.samline!=null ? r1.samline : new SamLine(r1, 0));
					if(r2!=null) {r2.samline=(r2.samline!=null ? r2.samline : new SamLine(r2, 1));}
				}
				list.add(r1);
			}
		}
		cros.add(list, num);
		return ro;
	}

	/** Expects a list of mixed mapped and unmapped reads; 
	 * only writes unmapped reads with unmapped mates */
	private int writeUnmappedOnly(ArrayList<Read> reads, long num) {
		if(reads==null || crosu==null || (reads.isEmpty() && !ordered)) {
			assert(!ordered || reads!=null);
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
		crosu.add(list, num);
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
	private final ConcurrentReadOutputStream cros;
	private final ConcurrentReadOutputStream crosu;
	
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
