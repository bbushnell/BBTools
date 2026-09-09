package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import aligner.SingleStateAlignerFlat2;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Orientation gate for reference-family sequence corpora: aligns every input sequence
 * against a reference (consensus) library in BOTH orientations and writes an oriented
 * corpus in which every retained record matches the library in the FORWARD direction.
 *
 * Motivation (rrnafallback STATUS cont.13-15, 2026-09-05): a handful of RefSeq GFF rRNA
 * rows carry inverted strand fields, so faithful extraction yields reverse-complemented
 * genes; three such records silently formed a reversed consensus model whose caller hits
 * landed on the wrong strand genome-wide. This tool is the flip-with-provenance
 * production design that prevents that class: reversed records are CORRECTED and tagged,
 * ambiguous records are EXCLUDED loudly, and nothing is ever silently kept.
 *
 * Per-record outcomes:
 *   FORWARD  - forward identity >= minid:            written unchanged.
 *   FLIPPED  - forward < minid but revcomp >= minid: written reverse-complemented, with
 *              " orientation=flipped fwdid=X rcid=Y" appended to the header (provenance).
 *   EXCLUDED - both orientations < minid:            written to outbad= (never kept).
 * Conservation is enforced with hard checks (not asserts): forward+flipped+excluded must
 * equal the input count, and the outcome ledger (out.tsv-style, one row per input) must
 * conserve input IDs.
 *
 * Uses SingleStateAlignerFlat2 (forward-only), the same aligner as MergeRibo's
 * acceptance gate, so minid here is directly comparable to MergeRibo minid.
 *
 * @author G11
 * @date September 5, 2026
 */
public class OrientSequences {

	public static void main(String[] args){
		Timer t=new Timer();
		OrientSequences x=new OrientSequences(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public OrientSequences(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in")){
				in=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("outbad")){
				outbad=b;
			}else if(a.equals("outledger") || a.equals("ledger")){
				outledger=b;
			}else if(a.equals("minid")){
				minID=Float.parseFloat(b);
			}else if(a.equals("ow") || a.equals("overwrite")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.setThreads(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if(in==null || ref==null || out==null){
			throw new RuntimeException("OrientSequences requires in=, ref=, and out=.");
		}
	}

	void process(Timer t){
		final ArrayList<Read> queries=loadReads(in);
		final ArrayList<Read> refs=loadReads(ref);
		if(queries.isEmpty()){throw new RuntimeException("No input sequences loaded from "+in);}
		if(refs.isEmpty()){throw new RuntimeException("No reference sequences loaded from "+ref);}
		outstream.println("Loaded "+queries.size()+" inputs, "+refs.size()+" references.");

		//One outcome slot per input, addressed by input index; threads write disjoint
		//slots so output preserves input order with no synchronization.
		final int n=queries.size();
		outcomes=new byte[n];
		fwdID=new float[n];
		rcID=new float[n];
		bestRefName=new String[n];
		final AtomicInteger next=new AtomicInteger(0);
		final int threads=Tools.mid(1, Shared.threads(), n);
		Thread[] workers=new Thread[threads];
		for(int i=0; i<threads; i++){
			workers[i]=new Thread(){
				@Override
				public void run(){
					final SingleStateAlignerFlat2 ssa=new SingleStateAlignerFlat2();
					for(int q=next.getAndIncrement(); q<n; q=next.getAndIncrement()){
						classify(queries.get(q), refs, ssa, q);
					}
				}
			};
			workers[i].start();
		}
		for(Thread w : workers){
			while(w.isAlive()){
				try{w.join();}catch(InterruptedException e){}
			}
		}

		//Write the three output streams plus the ledger, in input order.
		ByteStreamWriter bswOut=new ByteStreamWriter(FileFormat.testOutput(out, FileFormat.FASTA, null, true, overwrite, false, false));
		bswOut.start();
		ByteStreamWriter bswBad=(outbad==null ? null :
			new ByteStreamWriter(FileFormat.testOutput(outbad, FileFormat.FASTA, null, true, overwrite, false, false)));
		if(bswBad!=null){bswBad.start();}
		ByteStreamWriter bswLedger=(outledger==null ? null :
			new ByteStreamWriter(FileFormat.testOutput(outledger, FileFormat.TXT, null, true, overwrite, false, false)));
		if(bswLedger!=null){
			bswLedger.start();
			bswLedger.println("#qid\tqlen\tfwd_id\trc_id\tbest_ref\toutcome\tminid="+minID);
		}

		long nFwd=0, nFlip=0, nBad=0;
		long outRecords=0, badRecords=0;
		final ByteBuilder bb=new ByteBuilder();
		for(int q=0; q<n; q++){
			final Read r=queries.get(q);
			final byte outcome=outcomes[q];
			if(outcome==FORWARD){
				nFwd++; outRecords++;
				writeFasta(bswOut, r.id, r.bases, bb);
			}else if(outcome==FLIPPED){
				nFlip++; outRecords++;
				final byte[] rc=AminoAcid.reverseComplementBases(r.bases);
				//Provenance travels IN the header so downstream consumers can never mistake
				//a corrected record for an original one.
				writeFasta(bswOut, r.id+" orientation=flipped fwdid="+Tools.format("%.4f", fwdID[q])+
					" rcid="+Tools.format("%.4f", rcID[q]), rc, bb);
			}else{
				nBad++;
				if(bswBad!=null){badRecords++; writeFasta(bswBad, r.id, r.bases, bb);}
			}
			if(bswLedger!=null){
				bb.clear();
				bb.append(r.id).tab().append(r.length()).tab();
				bb.append(fwdID[q], 4).tab().append(rcID[q], 4).tab();
				bb.append(bestRefName[q]).tab();
				bb.append(outcome==FORWARD ? "FORWARD" : outcome==FLIPPED ? "FLIPPED" : "EXCLUDED").nl();
				bswLedger.print(bb.toBytes());
			}
		}
		bswOut.poisonAndWait();
		if(bswBad!=null){bswBad.poisonAndWait();}
		if(bswLedger!=null){bswLedger.poisonAndWait();}

		//Hard conservation checks (crash loud, never silently drop; hold without -ea):
		//this tool exists because 3 silently-reversed records poisoned a model
		//(rrnafallback STATUS cont.13); its own bookkeeping must never repeat that class.
		//(1) outcome partition covers every input; (2) write-branch counters (logical
		//bookkeeping incremented beside each write call, NOT an independent re-read of the
		//output files) must match the partition. Per-ID integrity holds structurally (one
		//slot per input index, written once in input order); INDEPENDENT verification is
		//external, via the ledger and the committed smoke script's file-level assertions.
		if(nFwd+nFlip+nBad!=n){
			throw new RuntimeException("CONSERVATION FAILED: forward="+nFwd+" + flipped="+nFlip+
				" + excluded="+nBad+" != inputs="+n);
		}
		if(outRecords!=nFwd+nFlip || (bswBad!=null && badRecords!=nBad)){
			throw new RuntimeException("WRITE-COUNT CONSERVATION FAILED: wrote "+outRecords+
				" oriented + "+badRecords+" excluded records for partition "+nFwd+"/"+nFlip+"/"+nBad);
		}
		if(nBad>0 && outbad==null){
			outstream.println("WARNING: "+nBad+" sequences failed both orientations and were EXCLUDED; "+
				"pass outbad= to capture them (they are never silently kept).");
		}

		t.stop();
		outstream.println("Inputs:  \t"+n);
		outstream.println("Forward: \t"+nFwd);
		outstream.println("Flipped: \t"+nFlip);
		outstream.println("Excluded:\t"+nBad);
		outstream.println("Time:    \t"+t);
	}

	/** Classify one query: forward-vs-library and revcomp-vs-library best identities. */
	private void classify(Read r, ArrayList<Read> refs, SingleStateAlignerFlat2 ssa, int index){
		float fwd=-1, rc=-1;
		String bestRef="-";
		final byte[] rcBases=AminoAcid.reverseComplementBases(r.bases);
		for(Read c : refs){
			final float f=ssa.align(r.bases, c.bases), v=ssa.align(rcBases, c.bases);
			if(Tools.max(f, v)>Tools.max(fwd, rc)){bestRef=c.id;}
			fwd=Tools.max(fwd, f);
			rc=Tools.max(rc, v);
		}
		bestRefName[index]=bestRef;
		fwdID[index]=fwd;
		rcID[index]=rc;
		//Forward wins ties deliberately: an original record is only modified when the
		//forward orientation FAILS the gate outright and revcomp passes it.
		outcomes[index]=(fwd>=minID ? FORWARD : rc>=minID ? FLIPPED : EXCLUDED);
	}

	private static void writeFasta(ByteStreamWriter bsw, String id, byte[] bases, ByteBuilder bb){
		synchronized(bsw){
			bb.clear();
			bb.append('>').append(id).nl();
			//Wrap at 70 columns to match the surrounding corpus files.
			for(int i=0; i<bases.length; i+=70){
				bb.append(bases, i, Tools.min(70, bases.length-i)).nl();
			}
			bsw.print(bb.toBytes());
		}
	}

	private ArrayList<Read> loadReads(String path){
		FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, false, false);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start();
		ArrayList<Read> list=new ArrayList<Read>();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){if(r.bases!=null){list.add(r);}}
			cris.returnList(ln);
		}
		ReadWrite.closeStream(cris);
		return list;
	}

	private static final byte FORWARD=0, FLIPPED=1, EXCLUDED=2;

	private String in=null;
	private String ref=null;
	private String out=null;
	private String outbad=null;
	private String outledger=null;
	private float minID=0.62f;//MergeRibo's default acceptance threshold, for direct comparability
	private boolean overwrite=true;
	private boolean verbose=false;

	private byte[] outcomes;
	private float[] fwdID;
	private float[] rcID;
	private String[] bestRefName;

	private PrintStream outstream=System.err;
}
