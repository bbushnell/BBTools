package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import aligner.SingleStateAlignerFlat2;
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
import tax.GiToTaxid;

/**
 * Per-sequence alignment acceptance diagnostic: aligns EVERY input sequence
 * against every record of a reference (consensus) library and emits exactly one
 * TSV row per input sequence with its best-matching reference and identity.
 *
 * Uses SingleStateAlignerFlat2 — the same aligner MergeRibo uses for its minID
 * acceptance gate — so identities here are directly comparable to MergeRibo's
 * minid threshold. Unlike MergeRibo, this tool NEVER groups by taxid and NEVER
 * reduces: taxid is a metadata column only, and the output row count must equal
 * the input sequence count. This measures per-sequence alignment acceptance,
 * NOT caller recall (caller recall = CallGenes output graded against genome-level
 * truth, a separate harness).
 *
 * Optional minlen/maxlen/maxns flags are recorded as a filter-verdict column
 * (what a MergeRibo-style pre-consensus filter WOULD do) but are never applied;
 * every sequence is aligned and reported regardless.
 *
 * @author G11
 * @date September 5, 2026
 */
public class AlignToConsensus {

	public static void main(String[] args){
		Timer t=new Timer();
		AlignToConsensus x=new AlignToConsensus(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public AlignToConsensus(String[] args){
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
			}else if(a.equals("minid")){
				minID=Float.parseFloat(b);
			}else if(a.equals("rcomp") || a.equals("rc")){
				reportRcomp=Parse.parseBoolean(b);
			}else if(a.equals("minlen")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("maxlen")){
				maxlen=Integer.parseInt(b);
			}else if(a.equals("maxns")){
				maxns=Integer.parseInt(b);
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
			throw new RuntimeException("AlignToConsensus requires in=, ref=, and out=.");
		}
	}

	void process(Timer t){
		final ArrayList<Read> queries=loadReads(in);
		final ArrayList<Read> refs=loadReads(ref);
		if(queries.isEmpty()){throw new RuntimeException("No query sequences loaded from "+in);}
		if(refs.isEmpty()){throw new RuntimeException("No reference sequences loaded from "+ref);}
		outstream.println("Loaded "+queries.size()+" queries, "+refs.size()+" references.");

		//One result slot per input, addressed by input index; threads write disjoint
		//slots so output preserves input order and no synchronization is needed.
		final byte[][] rows=new byte[queries.size()][];
		passFlags=new boolean[queries.size()];
		final AtomicInteger next=new AtomicInteger(0);
		final int threads=Tools.mid(1, Shared.threads(), queries.size());
		Thread[] workers=new Thread[threads];
		for(int i=0; i<threads; i++){
			workers[i]=new Thread(){
				@Override
				public void run(){
					final SingleStateAlignerFlat2 ssa=new SingleStateAlignerFlat2();
					final ByteBuilder bb=new ByteBuilder();
					for(int q=next.getAndIncrement(); q<queries.size(); q=next.getAndIncrement()){
						rows[q]=makeRow(queries.get(q), refs, ssa, bb, q);
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

		//Ledger completeness: exactly one row per input sequence, no reduction, and each
		//row carries its OWN query's id — count parity alone can't catch a duplicated or
		//swapped row (this tool exists because a taxid-grouped proxy silently reduced
		//8 inputs to 4 rows; see rrnafallback STATUS.md cont.6 retraction). These are
		//hard checks, not asserts: the ledger guarantee must hold even without -ea.
		for(int q=0; q<rows.length; q++){
			final String qid=queries.get(q).id;
			if(rows[q]==null){
				throw new RuntimeException("Missing result row for query index "+q+" ("+qid+
					"); every input must produce exactly one row.");
			}
			//Require the id to be terminated by the column tab so a row whose id merely
			//shares a PREFIX with the expected id (e.g. seq vs seq2) cannot pass.
			if(!Tools.startsWith(rows[q], qid) || rows[q].length<=qid.length() || rows[q][qid.length()]!='\t'){
				throw new RuntimeException("Row "+q+" does not begin with its exact query id '"+qid+
					"' followed by a tab; input-ID conservation violated.");
			}
		}

		ByteStreamWriter bsw=new ByteStreamWriter(FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, false, false));
		bsw.start();
		bsw.println("#qid\ttaxid\tqlen\tns\tfilter\tbest_ref\tbest_id\tverdict"+
			(reportRcomp ? "\trc_best_id\torient_flag" : "")+"\tminid="+minID);
		long pass=0;
		for(byte[] row : rows){bsw.print(row);}
		for(int q=0; q<rows.length; q++){if(passFlags[q]){pass++;}}
		bsw.poisonAndWait();

		t.stop();
		outstream.println("Queries:\t"+queries.size());
		outstream.println("Pass:\t"+pass+" ("+Tools.format("%.2f", pass*100.0/queries.size())+"%)");
		outstream.println("Fail:\t"+(queries.size()-pass));
		outstream.println("Time:\t"+t);
	}

	/** Build one TSV row for a query: best identity over ALL references, plus metadata. */
	private byte[] makeRow(Read r, ArrayList<Read> refs, SingleStateAlignerFlat2 ssa, ByteBuilder bb, int index){
		float bestID=-1;
		String bestRef="-";
		for(Read c : refs){
			float id=ssa.align(r.bases, c.bases);
			if(id>bestID){bestID=id; bestRef=c.id;}
		}
		assert(bestID>=0 && bestID<=1) : "Identity "+bestID+" outside [0,1] for "+r.id+
			" vs "+bestRef+"; SingleStateAlignerFlat2.align returns a fraction.";

		//Read.countNocalls counts exact 'N' bytes — the SAME semantics MergeRibo's maxns
		//filter uses (MergeRibo: r.countNocalls()>maxns), so the verdict column mirrors it exactly.
		final int ns=r.countNocalls();
		final String filter=(r.length()<minlen ? "short" : r.length()>maxlen ? "long" : ns>maxns ? "ns" : "-");
		final Integer tid=GiToTaxid.parseTaxidNumber(r.id, '|');
		final boolean passed=(bestID>=minID);
		passFlags[index]=passed;//disjoint index per thread; no synchronization needed

		bb.clear();
		bb.append(r.id).tab();
		bb.append(tid==null || tid.intValue()<0 ? "-" : tid.toString()).tab();//GiToTaxid returns -1, not null, when absent
		bb.append(r.length()).tab();
		bb.append(ns).tab();
		bb.append(filter).tab();
		bb.append(bestRef).tab();
		bb.append(bestID, 4).tab();
		bb.append(passed ? "PASS" : "FAIL");
		if(reportRcomp){
			//Orientation diagnostic (appended so existing column consumers are unaffected):
			//a FAIL whose revcomp passes is the reversed-record signature that caused the
			//cont.13 strand-flip incident; flag it REVERSED? so it can never hide again.
			float rcBest=-1;
			final byte[] rcBases=dna.AminoAcid.reverseComplementBases(r.bases);
			for(Read c : refs){rcBest=Tools.max(rcBest, ssa.align(rcBases, c.bases));}
			bb.tab().append(rcBest, 4).tab();
			bb.append(!passed && rcBest>=minID ? "REVERSED?" : "-");
		}
		bb.nl();
		return bb.toBytes();
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

	private String in=null;
	private String ref=null;
	private String out=null;
	private float minID=0.62f;//MergeRibo's default acceptance threshold, for direct comparability
	private boolean reportRcomp=false;//rcomp=t appends rc_best_id + orient_flag columns (REVERSED? on rc-passing FAILs)
	private int minlen=0;
	private int maxlen=Integer.MAX_VALUE;
	private int maxns=Integer.MAX_VALUE;
	private boolean overwrite=true;
	private boolean verbose=false;
	private boolean[] passFlags;

	private PrintStream outstream=System.err;
}
