package template;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.Writer;
import structures.ListNum;

/** Finite direct-worker validation checks using valid nucleotide reads.
 * Tests validation before writer handoff, not scheduling or malformed-input recovery.
 * @author Shinobu
 */
public final class TemplateStreamerValidationCase{

	/** Runs paired/unpaired cases with deferred or already completed validation. */
	public static void main(final String[] args){
		check(args.length==1 && ("expect=baseline".equals(args[0]) ||
				"expect=fixed".equals(args[0])), "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		final boolean oldConstructor=Read.VALIDATE_IN_CONSTRUCTOR;
		final boolean oldUpper=Read.TO_UPPER_CASE;
		final boolean oldLower=Read.LOWER_CASE_TO_N;
		final boolean oldSkip=Read.SKIP_SLOW_VALIDATION;
		try{
			Read.VALIDATE_IN_CONSTRUCTOR=false;
			Read.TO_UPPER_CASE=true;
			Read.LOWER_CASE_TO_N=false;
			Read.SKIP_SLOW_VALIDATION=false;
			for(boolean paired : new boolean[]{false, true}){
				for(boolean validated : new boolean[]{false, true}){
					runCase(paired, validated, fixed);
				}
			}
		}finally{
			Read.VALIDATE_IN_CONSTRUCTOR=oldConstructor;
			Read.TO_UPPER_CASE=oldUpper;
			Read.LOWER_CASE_TO_N=oldLower;
			Read.SKIP_SLOW_VALIDATION=oldSkip;
		}
	}

	/** Checks literal counters, retained batch ID, validation flags and bases at handoff. */
	private static void runCase(final boolean paired, final boolean validated,
			final boolean fixed){
		final String name=(paired ? "paired" : "single")+"-"+
				(validated ? "prevalidated" : "deferred");
		final Read root=new Read("acgt".getBytes(StandardCharsets.US_ASCII), null, "root", 0, false);
		final Read mate=paired ? new Read("tgca".getBytes(StandardCharsets.US_ASCII),
				null, "mate", 0, false) : null;
		if(mate!=null){root.mate=mate; mate.mate=root; mate.setPairnum(1);}
		if(validated){
			root.validate(true);
			if(mate!=null){mate.validate(true);}
		}
		final BatchSource source=new BatchSource(root);
		final RecordingWriter writer=new RecordingWriter();
		final A_SampleStreamerMT.ProcessThread worker=
				new A_SampleStreamerMT.ProcessThread(source, writer, 0);
		worker.run();
		final int reads=paired ? 2 : 1;
		check(worker.success && source.calls==2, name+" did not consume one batch and terminal null");
		check(worker.readsProcessedT==reads && worker.readsOutT==reads,
				name+" read counts: "+worker.readsProcessedT+","+worker.readsOutT);
		check(worker.basesProcessedT==4*reads && worker.basesOutT==4*reads,
				name+" base counts: "+worker.basesProcessedT+","+worker.basesOutT);
		check(writer.calls==1 && writer.id==0 && writer.root==root,
				name+" changed the batch, root or handoff count");
		final boolean expectedValidated=fixed || validated;
		check(writer.rootValidated==expectedValidated, name+" root validation at handoff");
		check(writer.rootBases.equals(expectedValidated ? "ACGT" : "acgt"),
				name+" root bases at handoff: "+writer.rootBases);
		if(paired){
			check(writer.mate==mate && writer.mateValidated==expectedValidated,
					name+" mate identity or validation at handoff");
			check(writer.mateBases.equals(expectedValidated ? "TGCA" : "tgca"),
					name+" mate bases at handoff: "+writer.mateBases);
		}else{check(writer.mate==null, name+" unexpectedly acquired a mate");}
		System.out.println("PASS\t"+name+"\t"+writer.rootValidated+"\t"+
				(paired ? Boolean.toString(writer.mateValidated) : "NA")+"\t"+
				writer.rootBases+"\t"+(paired ? writer.mateBases : "NA")+"\t"+reads+"\t"+(4*reads));
	}

	/** Fails with a concrete fixture observation independently of assertion settings. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Supplies one root batch, then terminal null. */
	private static final class BatchSource implements Streamer{
		BatchSource(final Read root_){root=root_;}
		@Override public ListNum<Read> nextList(){
			if(calls++>0){return null;}
			final ArrayList<Read> reads=new ArrayList<Read>();
			reads.add(root);
			return new ListNum<Read>(reads, 0);
		}
		@Override public String fname(){return "valid-literal-reads";}
		@Override public void start(){}
		@Override public void close(){}
		@Override public boolean paired(){return root.mate!=null;}
		@Override public int pairnum(){return 0;}
		@Override public long readsProcessed(){return root.pairCount();}
		@Override public long basesProcessed(){return root.pairLength();}
		@Override public void setSampleRate(final float rate, final long seed){throw new UnsupportedOperationException();}
		@Override public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException();}
		@Override public boolean hasMore(){return calls==0;}
		@Override public boolean errorState(){return false;}
		final Read root;
		int calls;
	}

	/** Captures flags and immutable base strings at submission, before later changes. */
	private static final class RecordingWriter implements Writer{
		@Override public void add(final ArrayList<Read> reads, final long id_){
			check(reads.size()==1, "Expected one submitted root");
			calls++;
			id=id_;
			root=reads.get(0);
			mate=root.mate;
			rootValidated=root.validated();
			rootBases=new String(root.bases, StandardCharsets.US_ASCII);
			if(mate!=null){
				mateValidated=mate.validated();
				mateBases=new String(mate.bases, StandardCharsets.US_ASCII);
			}
		}
		@Override public void addReads(final ListNum<Read> reads){add(reads.list, reads.id);}
		@Override public void addLines(final ListNum<SamLine> lines){throw new UnsupportedOperationException();}
		@Override public void start(){}
		@Override public long readsWritten(){return root==null ? 0 : root.pairCount();}
		@Override public long basesWritten(){return root==null ? 0 : root.pairLength();}
		@Override public void poison(){}
		@Override public boolean waitForFinish(){return false;}
		@Override public boolean poisonAndWait(){return false;}
		@Override public void finishError(){throw new UnsupportedOperationException();}
		@Override public String fname(){return "handoff-observations";}
		@Override public boolean errorState(){return false;}
		@Override public boolean finishedSuccessfully(){return true;}
		Read root, mate;
		String rootBases, mateBases;
		boolean rootValidated, mateValidated;
		long id=-1;
		int calls;
	}
}
