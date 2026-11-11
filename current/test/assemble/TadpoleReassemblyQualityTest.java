package test.assemble;

import java.util.Arrays;
import java.util.BitSet;

import assemble.ErrorTracker;
import assemble.Tadpole;
import assemble.Tadpole1;
import assemble.Tadpole2;
import kmer.AbstractKmerTableSet;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;
import ukmer.Kmer;

/**
 * Checks the quality-coordinate contract at both concrete reassembly entry points.
 * The shared pass supplies forward and reverse-complemented sequence copies. Tiny
 * subclasses observe each call, test the actual quality-dependent threshold helper,
 * then delegate to the real concrete inner method. Zero depths prevent proposals
 * and table queries, isolating orientation from other correction heuristics.
 * Optional reverse-call exception injection checks quality restoration on failure.
 * This is not a real-read correction-output or accuracy test. Each case runs in
 * its own JVM because assembler construction changes shared configuration.
 * @author Brian Bushnell, Sucrose
 */
public final class TadpoleReassemblyQualityTest{

	/** Static fixture entry point only. */
	private TadpoleReassemblyQualityTest(){}

	/** Accepts in=existing.fa, k=31|63, quality=asymmetric|symmetric|none and throwreverse=t|f.
	 * Exception injection requires asymmetric qualities and defaults to false.
	 * The none option selects a null array and retains null as the result label;
	 * PreParser rewrites a literal null option value before this parser sees it. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, TadpoleReassemblyQualityTest.class, false);
		try{
			String input=null, quality=null;
			int k=0;
			boolean throwReverse=false;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected named fixture options: "+arg);}
				if(pair[0].equalsIgnoreCase("in")){input=pair[1];}else if(pair[0].equalsIgnoreCase("k")){k=Integer.parseInt(pair[1]);}else if(pair[0].equalsIgnoreCase("quality")){quality="none".equals(pair[1]) ? "null" : pair[1];}else if(pair[0].equalsIgnoreCase("throwreverse")){throwReverse=Parse.parseBoolean(pair[1]);}else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			if(input==null || input.length()==0 || (k!=31 && k!=63) ||
					!("asymmetric".equals(quality) || "symmetric".equals(quality) || "null".equals(quality))){
				throw new IllegalArgumentException("Required in=existing.fa k=31|63 quality=asymmetric|symmetric|none.");
			}
			if(throwReverse && !"asymmetric".equals(quality)){
				throw new IllegalArgumentException("throwreverse=t requires quality=asymmetric to expose failed restoration.");
			}
			run(input, k, quality, throwReverse);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Executes one shared pass without loading reads, starting workers or opening outputs. */
	private static void run(final String input, final int k, final String quality, final boolean throwReverse){
		assert(k==31 || k==63) : "The fixture selects the two concrete kmer representations.";
		final int length=2*k+9;
		final byte[] bases=new byte[length];
		Arrays.fill(bases, (byte)'A');
		bases[length-1]='C';
		final byte[] quals="null".equals(quality) ? null : new byte[length];
		if(quals!=null){
			Arrays.fill(quals, (byte)20);
			if("asymmetric".equals(quality)){quals[k+2]=0; quals[k+6]=40;}
		}
		final Trace trace=new Trace(k, quality, bases, quals, throwReverse);
		Kmer.PACKED=false;
		AbstractKmerTableSet.MASK_CORE=k<=31;
		Shared.COMMAND_LINE=null;
		final String[] options={"in="+input, "out=null", "k="+k, "mode=correct", "t=1",
				"prefilter=0", "prealloc=f", "ways=1", "initialsize=128", "minprob=0",
				"ecc=f", "ecco=f", "merge=f", "markerrors=f", "hpindel=f", "fixindels=f",
				"hollow=f", "hollowbranchstop=f", "hollowbranchrollback=f", "hollowrescue=f",
				"aggressive=f", "conservative=f", "deadzone=0", "errorpath=1",
				"errormult1=16", "errormult2=2.6", "errormultqfactor=0.002",
				"errorlowerconst=4", "mincountcorrect=3", "showstats=f"};
		final Tadpole tad=k==31 ? new ShortProbe(options, trace) : new LongProbe(options, trace);
		try{
			System.err.println("QUALITY_ORIENTATION_START k="+k+" quality="+quality+" length="+length+
					" probe_index="+trace.position+" assertions="+TadpoleReassemblyQualityTest.class.desiredAssertionStatus());
			final IntList counts=new IntList(length-k+1), counts2=new IntList(length-k+1);
			for(int i=0; i<length-k+1; i++){counts.add(0);}
			final ErrorTracker tracker=new ErrorTracker();
			int corrected=-1;
			boolean caught=false;
			try{
				corrected=tad.reassemble_pass(bases, quals, new ByteBuilder(length),
					new ByteBuilder(length), new int[4], counts, counts2, tracker, 1,
					new Kmer(k), new Kmer(k), new BitSet(length));
			}catch(RuntimeException e){
				if(e!=trace.injectedFailure){throw e;}
				caught=true;
			}
			check(caught==throwReverse, "Reverse-inner exception was swallowed or injected in an ordinary case.");
			check(trace.calls==2 && trace.completed==(throwReverse ? 1 : 2),
					"Unexpected observed/completed concrete inner calls.");
			check(corrected==(throwReverse ? -1 : 0) && tracker.detected()==0 && tracker.corrected()==0,
					"Zero-depth fixture unexpectedly proposed or installed a correction.");
			check(Arrays.equals(bases, trace.originalBases), "Shared pass changed original bases in the no-proposal control.");
			final boolean qualitiesPreserved=Arrays.equals(quals, trace.originalQuals);
			check(counts.size==length-k+1, "Shared pass changed the original count-list size.");
			for(int i=0; i<counts.size; i++){check(counts.get(i)==0, "Shared pass changed an original depth.");}
			final boolean passed=trace.aligned && trace.decisions && qualitiesPreserved;
			if(throwReverse){
				System.out.println("QUALITY_RESTORATION_RESULT k="+k+" quality="+quality+
						" pass="+passed+" calls=2 completed=1 exception_preserved=true bases_preserved=true qualities_preserved="+
						qualitiesPreserved+" depths_preserved=true");
				check(passed, "QUALITY_RESTORATION_CONTRACT_MISMATCH");
				return;
			}
			check(qualitiesPreserved, "Shared pass changed original qualities in the no-proposal control.");
			System.out.println("QUALITY_ORIENTATION_RESULT k="+k+" quality="+quality+
					" pass="+passed+" calls=2 corrections=0 bases_preserved=true qualities_preserved=true depths_preserved=true");
			check(passed, "QUALITY_ORIENTATION_CONTRACT_MISMATCH");
		}finally{tad.clearData();}
	}

	/** Always-on checks distinguish setup failures from the named orientation contract failure. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Records one process's two calls against independently constructed literal expectations. */
	private static final class Trace{
		Trace(final int k_, final String quality_, final byte[] bases, final byte[] quals, final boolean throwReverse){
			k=k_; quality=quality_; position=k+2;
			originalBases=bases.clone(); originalQuals=quals==null ? null : quals.clone();
			injectedFailure=throwReverse ? new RuntimeException("QUALITY_REVERSE_INNER_INJECTED_FAILURE") : null;
			check(originalBases.length-1-position==k+6, "Fixture quality positions must mirror one another.");
		}

		/** Captures arguments before the real inner method can change its sequence buffer. */
		void observe(final ByteBuilder bases, final byte[] quals, final IntList counts, final boolean decision){
			check(calls<2 && bases.length()==originalBases.length, "Unexpected inner-call count or sequence length.");
			final boolean reverse=calls==1;
			for(int i=0; i<bases.length(); i++){
				final byte expected=reverse ? (byte)(i==0 ? 'G' : 'T') : originalBases[i];
				check(bases.get(i)==expected, "Inner call received an unexpected sequence orientation.");
			}
			check(counts.size==bases.length()-k+1, "Inner call received an invalid profile size.");
			for(int i=0; i<counts.size; i++){check(counts.get(i)==0, "Inner call received a nonzero depth.");}
			check((quals==null)==(originalQuals==null) && (quals==null || quals.length==bases.length()),
					"Inner call changed quality presence or length.");
			boolean matches=true;
			if(quals!=null){for(int i=0; i<quals.length; i++){
				matches&=quals[i]==originalQuals[reverse ? quals.length-1-i : i];
			}}
			final int actualQ=quals==null ? 20 : quals[position];
			final int expectedQ="asymmetric".equals(quality) ? (reverse ? 40 : 0) : 20;
			final boolean expectedDecision="asymmetric".equals(quality) && !reverse;
			aligned&=matches; decisions&=decision==expectedDecision;
			System.out.println("QUALITY_ORIENTATION_CALL k="+k+" quality="+quality+
					" direction="+(reverse ? "reverse" : "forward")+" probe_index="+position+
					" original_index="+(reverse ? bases.length()-1-position : position)+" aligned="+matches+
					" observed_q="+actualQ+" expected_q="+expectedQ+" decision="+decision+
					" expected_decision="+expectedDecision);
			calls++;
			// Throw only after the reverse arguments have been observed, before concrete inner work.
			if(reverse && injectedFailure!=null){throw injectedFailure;}
		}

		final int k, position;
		final String quality;
		final byte[] originalBases, originalQuals;
		/** Exact object propagated through the shared pass; unrelated exceptions remain setup errors. */
		final RuntimeException injectedFailure;
		int calls, completed;
		boolean aligned=true, decisions=true;
	}

	/** Observes the shared caller's arguments and delegates to the short-word implementation. */
	private static final class ShortProbe extends Tadpole1{
		ShortProbe(final String[] args, final Trace trace_){super(args, false); trace=trace_;}
		@Override
		public int reassemble_inner(final ByteBuilder bases, final byte[] quals, final int[] rightCounts,
				final IntList counts, final int extension, final Kmer kmer, final Kmer regenKmer){
			final byte q=quals==null ? 20 : quals[trace.position];
			trace.observe(bases, quals, counts, isError(82, 5, q));
			final int result=super.reassemble_inner(bases, quals, rightCounts, counts, extension, kmer, regenKmer);
			check(result==0, "Short concrete inner unexpectedly corrected the zero-depth fixture.");
			trace.completed++;
			return result;
		}
		final Trace trace;
	}

	/** Observes the shared caller's arguments and delegates to the packed-word implementation. */
	private static final class LongProbe extends Tadpole2{
		LongProbe(final String[] args, final Trace trace_){super(args, false); trace=trace_;}
		@Override
		public int reassemble_inner(final ByteBuilder bases, final byte[] quals, final int[] rightCounts,
				final IntList counts, final int extension, final Kmer kmer, final Kmer regenKmer){
			final byte q=quals==null ? 20 : quals[trace.position];
			trace.observe(bases, quals, counts, isError(82, 5, q));
			final int result=super.reassemble_inner(bases, quals, rightCounts, counts, extension, kmer, regenKmer);
			check(result==0, "Long concrete inner unexpectedly corrected the zero-depth fixture.");
			trace.completed++;
			return result;
		}
		final Trace trace;
	}
}
