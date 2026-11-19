package test.cardinality;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import kmer.AbstractKmerTable;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData;
import prok.NcrnaTextIO.Output;
import prok.RrnaResourceIO;
import shared.KillSwitch;
import shared.Shared;
import structures.ByteBuilder;

/** Matched native versus primitive-consumer fill timing on saved sequences.
 * The actual forward masked/tagged dual fill is a separate workload, not a
 * claim of equivalent native canonical single-k semantics.
 * @author Ganyu
 */
public final class DdlFillTiming implements CardinalityTracker.KmerConsumer{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, DdlFillTiming.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String panels=null;
			String selection=null;
			int seconds=10;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("panels")){panels=b;}
				else if(a.equals("selection")){selection=b;}
				else if(a.equals("seconds")){seconds=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown fill-timing argument: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions && parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1 && Shared.threads()==1
				&& (seconds==0 || seconds>=10) && ((panels==null)==(selection==null)),
				"Use out=fresh.tsv t=1 -ea; seconds=0 finite fixture or >=10 with complete saved panels/selection");
			require(seconds==0 || panels!=null, "Measured timing requires the saved-window sequence population");
			require(!DynamicDemiLog.blacklistExists(), "Fill timing isolates insertion without an unrequested global blacklist");
			RrnaResourceIO.fresh(parser.out1, "");
			final int exponent=DynamicDemiLog.exponentBits();
			DynamicDemiLog.setExponent(4);
			try(Output output=new Output(parser.out1)){
				final DdlFillTiming test=new DdlFillTiming();
				test.verify(new byte[][]{
					new byte[0], "ACG".getBytes(StandardCharsets.US_ASCII),
					"ACGTACGTACGTACGT".getBytes(StandardCharsets.US_ASCII),
					"acgtaNNNacgtacgtacgtNNacgta".getBytes(StandardCharsets.US_ASCII),
					"AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT".getBytes(StandardCharsets.US_ASCII)});
				byte[][] sequences=null;
				if(panels!=null){
					final NcrnaDdlSavedData data=NcrnaDdlSavedData.load(panels, selection, 508, 3153);
					require(data.windowCount()==2531, "Timing must retain all 2531 original saved windows");
					sequences=new byte[data.windowCount()][];
					for(int i=0; i<sequences.length; i++){sequences[i]=data.windowBases(i);}
					test.verify(sequences);
				}
				output.writer.println("method\tround\twindows_per_sweep\tsequence_calls\tmain_thread_cpu_seconds\tcpu_microseconds_per_window\tkey_checksum\tsemantics\tscope");
				if(seconds>0){test.measure(sequences, seconds, output);}
				pp.outstream.println(seconds==0 ? "DDL_FILL_FINITE_PARITY_PASS" : "DDL_FILL_PARITY_AND_TIMING_PASS");
			}finally{DynamicDemiLog.setExponent(exponent);}
		}catch(Throwable failure){
			failure.printStackTrace(pp.outstream);
			KillSwitch.kill("DDL fill timing failed");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private DdlFillTiming(){
		require(clock.isCurrentThreadCpuTimeSupported(), "Timing requires the main thread's CPU clock");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
	}

	/** Complete key equality precedes every matched control timing. Native add()
	 * also maintains its added counter; only key filling is compared here. */
	private void verify(final byte[][] sequences){
		assert(sequences!=null && sequences.length>0) : "Parity uses a nonempty sequence roster";
		final DynamicDemiLog literal=DynamicDemiLog.create(512, 9, 12345, 0, false, false);
		for(final int code:new int[]{0,1}){
			final byte[] bases=new byte[9];
			Arrays.fill(bases, (byte)(code==0 ? 'A' : 'C'));
			fill(bases, 4);
			literal.clear();
			literal.hashAndStore((code*((1L<<14)-1)/3)|(1L<<14));
			literal.hashAndStore(((code*((1L<<18)-1)/3)&middleMask)|(1L<<18));
			require(Arrays.equals(literal.maxArray(), active.maxArray()),
				"Literal homopolymers must retain separate leading-tag k7 and middle-masked k9 insertions");
		}
		for(final byte[] bases:sequences){
			for(int pair=0; pair<2; pair++){
				fill(bases, pair*2);
				final char[] expected=active.maxArray().clone();
				fill(bases, pair*2+1);
				require(Arrays.equals(expected, active.maxArray()), "Native/shared iterator complete key mismatch at k="+(pair==0 ? 7 : 9));
			}
		}
		System.err.println("DDL_FILL_ALL_KEYS_PARITY_PASS\tsequences="+sequences.length+"\tk=7,9\tbins=512");
	}

	/** Reuses native objects and consumes every resulting key on both paths.
	 * Allocation, data loading and parity are outside the timed interval. */
	private void measure(final byte[][] sequences, final int seconds, final Output output){
		assert(sequences!=null && sequences.length==2531 && seconds>=10)
			: "Saved-window timing requires the full fixed population and at least ten CPU seconds";
		final long minimum=seconds*1_000_000_000L;
		final long[] expected=new long[NAMES.length];
		for(int arm=0; arm<NAMES.length; arm++){
			for(final byte[] bases:sequences){fill(bases, arm); expected[arm]+=checksum();}
		}
		require(expected[0]==expected[1] && expected[2]==expected[3], "Matched native/iterator controls must have equal whole-sweep checksums");
		final ByteBuilder row=new ByteBuilder();
		for(int round=0; round<3; round++){
			for(int order=0; order<NAMES.length; order++){
				final int arm=order==4 ? 4 : order^(round&1);
				long cpu=0;
				long calls=0;
				long sum=0;
				do{
					final long start=clock.getCurrentThreadCpuTime();
					for(final byte[] bases:sequences){fill(bases, arm); sum+=checksum();}
					cpu+=clock.getCurrentThreadCpuTime()-start;
					calls+=sequences.length;
				}while(cpu<minimum);
				require(sum==Math.multiplyExact(expected[arm], calls/sequences.length), "Timed fills must preserve the consumed complete key population");
				row.clear().append(NAMES[arm]).tab().append(round+1).tab().append(sequences.length).tab().append(calls)
					.tab().append(cpu/1e9, 9).tab().append(cpu/(1000.0*calls), 9).tab().append(sum)
					.tab().append(arm==4 ? "FORWARD_STRICT_K7_MASKED9_LEADING_BIT_MIXED512" : "CANONICAL_NATIVE_AMBIGUITY_SINGLE_K_512")
					.append("\tCLEAR_FILL_AND_ALL_KEY_CHECKSUM_IO_CONSTRUCTION_EXCLUDED\n");
				output.writer.print(row);
			}
		}
	}

	private void fill(final byte[] bases, final int arm){
		assert(bases!=null && arm>=0 && arm<NAMES.length) : "Each timed arm borrows one complete saved sequence";
		mode=arm;
		active=arm<2 ? shortDdl : longDdl;
		active.clear();
		if(arm==0 || arm==2){active.hash(bases, null);}
		else if(arm==1 || arm==3){
			final int k=arm==1 ? 7 : 9;
			CardinalityTracker.forEachKmer(bases, k, k, true, false, this);
		}else{
			// Match the word processing in the producer's two forward passes.
			// This sink does not measure its frequency/blacklist callback branches.
			CardinalityTracker.forEachKmer(bases, 7, 7, false, true, this);
			CardinalityTracker.forEachKmer(bases, 9, 9, false, true, this);
		}
	}

	@Override
	public void accept(final long word, final int k){
		assert(active!=null && (k==7 || k==9)) : "The consumer belongs to the active k7/k9 fill";
		active.hashAndStore(mode==4 ? ((word&(k==7 ? -1L : middleMask))|(1L<<(2*k))) : word);
	}

	private long checksum(){
		assert(active!=null) : "Consume the same filled object produced by the selected arm";
		long sum=0;
		for(final char value:active.maxArray()){sum+=value;}
		return sum;
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final String[] NAMES={"NATIVE_K7", "ITERATOR_HASHANDSTORE_K7", "NATIVE_K9", "ITERATOR_HASHANDSTORE_K9", "TAGGED_DUAL_REFERENCE_TWO_PASSES"};
	private final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
	private final DynamicDemiLog shortDdl=DynamicDemiLog.create(512, 7, 12345, 0, false, false);
	private final DynamicDemiLog longDdl=DynamicDemiLog.create(512, 9, 12345, 0, false, false);
	private final long middleMask=AbstractKmerTable.makeMiddleMask(9, false, 1);
	private DynamicDemiLog active;
	private int mode;
}
